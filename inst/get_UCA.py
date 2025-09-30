import argparse
import sys
from Bio.Seq import Seq
import json
import logomaker
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
import olga.load_model as load_model
import olga.generation_probability as pgen
import pandas as pd 
import random
import os

def read_starting_and_ending_points(file_path):
    with open(file_path, 'r') as file:
        content = file.read().strip()
        starting_point, ending_point = map(int, content.split())
    return starting_point, ending_point

def generate_codon_list(codons, starting_germline, starting_point, ending_point, cdr3_only, tree_df):
    if cdr3_only:
        codons_cdr3 = [(i, starting_germline[i:i+3]) for i in range(starting_point, ending_point, 3)]
        positions = [pos for pos, _ in codons_cdr3]
        positions = [
            pos for pos in positions
            if len(tree_df[tree_df['site'] == pos / 3]) > 1
        ]
        new_order = random.sample(positions, len(positions))
    else:
        new_order = random.sample(range(len(codons)), len(codons))
    return new_order

def propose_new_combination(current_germline, codon_number, igphyml_df):
    current_codon = current_germline[codon_number:codon_number+3]
    sub_df = igphyml_df[igphyml_df["site"] == codon_number/3]
    sub_df = sub_df[sub_df["codon"] != current_codon]
    new_codon = sub_df.sample(n=1)
    new_germline = current_germline[:codon_number] + new_codon["codon"].values[0] + current_germline[codon_number+3:]
    return new_germline 
    
def get_new_lhood(proposed_combination, starting_point, ending_point, pgen_model, igphyml_df):
    cdr3 = proposed_combination[starting_point:ending_point]
    new_pgen = pgen_model.compute_nt_CDR3_pgen(cdr3)
    if new_pgen == 0:
        new_pgen = new_pgen + 1e-323
    new_pgen = np.log(new_pgen)
    new_codons = [proposed_combination[i:i+3] for i in range(0, len(proposed_combination), 3)]
    partial_likelihoods = []
    for i in range(0, igphyml_df['site'].max()):
        selection = igphyml_df.loc[(igphyml_df['site'] == i) & (igphyml_df['codon'] == new_codons[i]), 'partial_likelihood']
        if selection.empty:
            matched_value = -1e6
        else:
            matched_value = selection.values[0]
        if not np.isfinite(matched_value):
            matched_value = -1e6
        partial_likelihoods.append(matched_value)
    new_tree_likelihood = sum(partial_likelihoods)
    new_lhood = new_pgen + new_tree_likelihood
    lhood_vector = [new_pgen, new_tree_likelihood, new_lhood]
    return lhood_vector

def translate_to_amino_acid(codon):
    return str(Seq(codon).translate())

def process_option(option, filtered_df, codon, current_germline, starting_point, ending_point, igphyml_df_new, pgen_model):
    testing_df = filtered_df.iloc[option]
    testing_germline = current_germline 
    testing_germline = [testing_germline[i:i+3] for i in range(0, len(testing_germline), 3)]
    testing_germline[codon // 3] = testing_df['codon']
    testing_germline = ''.join(testing_germline)
    new_lhoods = get_new_lhood(testing_germline, starting_point, ending_point, pgen_model, igphyml_df_new)
    temp_data = {
    'site': codon//3,
    'codon': testing_df['codon'],
    'joint_log_likelihood': new_lhoods[2],
    'tree_log_likelihood': new_lhoods[1],
    'pgen_log_likelihood': new_lhoods[0],  
    }
    return testing_germline, new_lhoods, temp_data

def process_option_wrapper(args):
    return process_option(*args)

def process_nt_option(base, codons, codon_numbers, codon_group_index, site_codon_indx, current_germline, starting_point, ending_point, igphyml_df_new, pgen_model):
    codon_list = list(codons[codon_group_index])
    codon_list[site_codon_indx] = base
    test_codon = ''.join(codon_list)
    if test_codon in ["TAA", "TAG", "TGA"]:
        return None, None, None
    seq_list = list(current_germline)
    for idx, b in zip(codon_numbers[codon_group_index], test_codon):
        seq_list[idx-1] = b
    testing_germline = ''.join(seq_list)
    new_lhoods = get_new_lhood(testing_germline, starting_point, ending_point, pgen_model, igphyml_df_new)
    temp_data = {
        'site': (max(codon_numbers[codon_group_index]) // 3) - 1,
        'nt_site': codon_numbers[codon_group_index][site_codon_indx],
        'codon': test_codon,
        'nt': base,
        'joint_log_likelihood': new_lhoods[2],
        'tree_log_likelihood': new_lhoods[1],
        'pgen_log_likelihood': new_lhoods[0],  
    }
    return testing_germline, new_lhoods, temp_data


def get_updated_germline(starting_germline, starting_cdr3, igphyml_df, pgen_model, starting_point, ending_point, cdr3_only=True, max_iter=100, nproc = 1, quiet = 0):
    current_codons = [starting_germline[i:i+3] for i in range(0, len(starting_germline), 3)]
    current_germline = starting_germline
    current_lhoods = get_new_lhood(current_germline, starting_point, ending_point, pgen_model, igphyml_df)
    starting_site = starting_point/3
    ending_site = ending_point/3 - 1
    starting_options = igphyml_df[igphyml_df["site"] == starting_site]
    starting_options = starting_options[starting_options['codon'].apply(translate_to_amino_acid) == 'C']
    ending_options = igphyml_df[igphyml_df["site"] == ending_site]
    ending_options = ending_options[ending_options['codon'].apply(translate_to_amino_acid).isin(['W', 'F'])]
    filtered_igphyml_df = igphyml_df[~igphyml_df['site'].isin([starting_site, ending_site])]
    igphyml_df_new = pd.concat([filtered_igphyml_df, starting_options, ending_options])
    igphyml_df_new = igphyml_df_new.sort_values(by='site', ascending=True)
    still_improving = True
    tested_combinations = []
    tested_lhoods_pgen = []
    tested_lhoods_tree = []
    tested_lhoods_combo = []
    tested_iterations = []
    iteration = 0
    
    while still_improving and iteration < max_iter:
        still_improving = False
        codon_list = generate_codon_list(current_codons, current_germline, starting_point, ending_point, cdr3_only, igphyml_df_new)
        iteration += 1
        iteration_data = []
        for codon in codon_list:           
            filtered_df = igphyml_df_new[igphyml_df_new["site"] == codon/3]
            with Pool(nproc) as pool:
                results = pool.map(process_option_wrapper, [(option, filtered_df, codon, current_germline, starting_point, ending_point, igphyml_df_new, pgen_model) for option in range(0, len(filtered_df))])
            results_array = np.array([result[1] for result in results])
            best_index = np.argmax(results_array[:, 2])
            tested_combinations.append(results[best_index][0])
            tested_lhoods_pgen.append(results[best_index][1][0])
            tested_lhoods_tree.append(results[best_index][1][1])
            tested_lhoods_combo.append(results[best_index][1][2])
            indices = [i for i, x in enumerate(codon_list) if x == codon]
            tested_iterations.extend([f"{iteration}_{index}" for index in indices])
            if results_array[best_index][2] > current_lhoods[2]:
                current_lhoods = list(results_array[best_index])
                current_germline = results[best_index][0]
                still_improving = True
            
            # Collect temp_data for this iteration
            for result in results:
                temp_data = result[2]
                iteration_data.append(temp_data)
        
        iteration_df = pd.DataFrame(iteration_data)
        iteration_df = iteration_df.sort_values(by=['site', 'codon'], ascending=[True, True])
        iteration_df['relative_likelihood'] = iteration_df.groupby('site')['joint_log_likelihood'].transform(
            lambda x: np.exp(x - np.logaddexp.reduce(x))
        )

         

    data = {
    'germline': tested_combinations,
    'pgen_lhood': tested_lhoods_pgen,
    'tree_lhood': tested_lhoods_tree,
    'combo_lhood': tested_lhoods_combo,
    'iteration': tested_iterations,
    }
    data = pd.DataFrame(data)
    return current_germline, current_lhoods, data, iteration_df

def get_updated_germline_nt(starting_germline, starting_cdr3, igphyml_df, pgen_model, starting_point, ending_point, cdr3_only=True, max_iter=100, nproc=1):
    current_codons = [starting_germline[i:i+3] for i in range(0, len(starting_germline), 3)]
    current_germline = starting_germline
    current_lhoods = get_new_lhood(current_germline, starting_point, ending_point, pgen_model, igphyml_df)
    starting_site = starting_point/3 + 1
    ending_site = ending_point/3 - 1
    starting_options = igphyml_df[igphyml_df["site"] == starting_site]
    starting_options = starting_options[starting_options['codon'].apply(translate_to_amino_acid) == 'C']
    ending_options = igphyml_df[igphyml_df["site"] == ending_site]
    ending_options = ending_options[ending_options['codon'].apply(translate_to_amino_acid).isin(['W', 'F'])]
    filtered_igphyml_df = igphyml_df[~igphyml_df['site'].isin([starting_site, ending_site])]
    igphyml_df_new = pd.concat([filtered_igphyml_df, starting_options, ending_options])
    igphyml_df_new = igphyml_df_new.sort_values(by='site', ascending=True)
    still_improving = True
    tested_combinations = []
    tested_lhoods_pgen = []
    tested_lhoods_tree = []
    tested_lhoods_combo = []
    tested_iterations = []
    iteration = 0
    junction_sites = list(range(starting_point + 1, ending_point + 1))
    junction_site_groups = [junction_sites[i:i+3] for i in range(0, len(junction_sites), 3)]
    positions = [max(group) // 3 - 1 for group in junction_site_groups]
    positions = [
        pos for pos in positions
        if len(igphyml_df_new[igphyml_df_new['site'] == pos]) > 1
    ]
    positions = [
        group for group in junction_site_groups
        if (max(group) // 3 - 1) in positions
    ]
    junction_sites = [site for group in positions for site in group]

    while still_improving and iteration < max_iter:
        still_improving = False
        current_codons = [current_germline[i:i+3] for i in range(0, len(current_germline), 3)]
        codon_numbers = [list(range(i, min(i+3, len(current_germline)+1))) for i in range(1, len(current_germline)+1, 3)]
        iteration += 1
        iteration_data = []
        random_junction_sites = random.sample(junction_sites, len(junction_sites))
        
        for i in random_junction_sites:
            codon_group_index = next(idx for idx, group in enumerate(codon_numbers) if i in group)
            site_codon_indx = codon_numbers[codon_group_index].index(i)
            testing_bases = ['A', 'C', 'G', 'T']
            results = []
            for base in testing_bases:
                result = process_nt_option(base, current_codons, codon_numbers, codon_group_index, site_codon_indx, current_germline, starting_point, ending_point, igphyml_df_new, pgen_model)
                if result is not None and result[0] is not None and result[1] is not None:  # Only append if not a stop codon
                    results.append(result)
            if not results:
                continue
            results_array = np.array([result[1] for result in results])
            if results_array.size == 0:
                continue
            best_index = np.argmax(results_array[:, 2])
            tested_combinations.append(results[best_index][0])
            tested_lhoods_pgen.append(results[best_index][1][0])
            tested_lhoods_tree.append(results[best_index][1][1])
            tested_lhoods_combo.append(results[best_index][1][2])
            indices = [i for i, x in enumerate(testing_bases) if x == base]
            tested_iterations.extend([f"{iteration}_{index}" for index in indices])
            if results_array[best_index][2] > current_lhoods[2]:
                current_lhoods = list(results_array[best_index])
                current_germline = results[best_index][0]
                still_improving = True

            # Collect temp_data for this iteration
            for result in results:
                temp_data = result[2]
                iteration_data.append(temp_data)

        iteration_df = pd.DataFrame(iteration_data)
        iteration_df = iteration_df.sort_values(by=['site', 'codon'], ascending=[True, True])
        iteration_df['relative_likelihood'] = iteration_df.groupby('site')['joint_log_likelihood'].transform(
           lambda x: np.exp(x - np.logaddexp.reduce(x))
        ) 

    data = {
    'germline': tested_combinations,
    'pgen_lhood': tested_lhoods_pgen,
    'tree_lhood': tested_lhoods_tree,
    'combo_lhood': tested_lhoods_combo,
    'iteration': tested_iterations,
    }
    data = pd.DataFrame(data)
    return current_germline, current_lhoods, data, iteration_df

def process_row(row):
    directory = row['directory']
    id = row['id']
    quiet = row['quiet']
    model_folder = row['model_folder']
    model_folder_igk = row['model_folder_igk']
    model_folder_igl = row['model_folder_igl']
    max_iters = row['max_iters']
    
    clone_number = row['clone_ids']
    base_string = directory + "/" + id + "_" + clone_number
    if quiet > 0:
        print("\nRunning on clone", clone_number, "with locus", row['chains'], "saved at", base_string)
    
    if not os.path.exists(base_string):
        os.makedirs(base_string)
    chain = row['chains']

    table_path = row['tree_table']    
    igphyml_df = pd.read_csv(table_path, sep = "\t", header = None)
    igphyml_df.columns = ["site", "codon", "partial_likelihood", "log_likelihood_site", "upper_partial_log_likelihood", "upper_partial_likelihood", "equilibrium"]
    igphyml_df['value'] = igphyml_df['partial_likelihood'] + np.log(igphyml_df['equilibrium'])
    germline_string = row['starting_germline']
    junction_string = row['junction_locations']
    starting_point, ending_point = read_starting_and_ending_points(junction_string)
    with open(germline_string, "r") as f:
        starting_germline = f.read().strip()
    starting_cdr3 = starting_germline[starting_point:ending_point]

    v_gene = starting_germline[0:starting_point]
    j_gene = starting_germline[ending_point:len(starting_germline)]
    if "N" in v_gene or "N" in j_gene:
        if(quiet > 0):
            print("N in the J gene -- adding to the codon df")
        codons = [starting_germline[i:i+3] for i in range(0, len(starting_germline), 3)]
        indices_with_N = [index for index, value in enumerate(codons) if 'N' in value and not (len(v_gene) // 3 + 1 <= index <= len(codons) - len(j_gene) // 3  - 1)]
        new_rows = []
        for index in indices_with_N:
            sub = igphyml_df[igphyml_df['site'] == index]
            sub_sum = sub['value'].sum()
            new_rows.append({'site': index, 'codon': codons[index], 'partial_likelihood': sub_sum, 'nope': 0, 'nada': 0, 'no': 0, 'equilibrium': 0, 'value': 0})
        igphyml_df = pd.concat([igphyml_df, pd.DataFrame(new_rows)], ignore_index=True)
        igphyml_df = igphyml_df.sort_values(by='site', ascending=True)
        cdr3 = starting_germline.split(v_gene)
        cdr3  = cdr3[1].split(j_gene)[0]
        starting_germline = v_gene + cdr3.replace("N", "C") + j_gene

    
    if "N" in starting_cdr3:
        starting_cdr3 = starting_cdr3.replace("N", "C")
        starting_germline = v_gene + starting_cdr3 + j_gene

    codons = [starting_germline[i:i+3] for i in range(0, len(starting_germline), 3)]

    if chain == "IGH":
        model_folder = model_folder
    elif chain == "IGK":
        model_folder = model_folder_igk
    elif chain == "IGL":
        model_folder = model_folder_igl
    else:
        print("Invalid chain type. Please use IGH, IGK, or IGL.")
        sys.exit(1)

    params_file_name = model_folder + '/model_params.txt'
    marginals_file_name = model_folder + '/model_marginals.txt'
    V_anchor_pos_file = model_folder + '/V_gene_CDR3_anchors.csv'
    J_anchor_pos_file = model_folder + '/J_gene_CDR3_anchors.csv'

    if chain == "IGH":
        genomic_data = load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
    else:
        genomic_data = load_model.GenomicDataVJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = load_model.GenerativeModelVJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = pgen.GenerationProbabilityVJ(generative_model, genomic_data)
    
    values = get_updated_germline_nt(starting_germline, starting_cdr3, igphyml_df, pgen_model, starting_point, ending_point, cdr3_only=True, max_iter=int(max_iters), nproc = 1)

    new_germline = values[0]
    new_lhoods = values[1]
    data = values[2]
    iteration_df = values[3]

    collapsed_df = iteration_df.drop_duplicates(subset=['site', 'codon']).copy()
    collapsed_df['amino_acid'] = collapsed_df['codon'].apply(translate_to_amino_acid)
    logo_df = collapsed_df.groupby(['site', 'amino_acid'], as_index=False).agg({
        'relative_likelihood': 'sum'
    })
    pwm = logo_df.pivot(index='site', columns='amino_acid', values='relative_likelihood').fillna(0)

    iteration_df['amino_acid'] = iteration_df['codon'].apply(translate_to_amino_acid)
    iteration_df['uca_codon'] = iteration_df.groupby('site')['relative_likelihood'].transform(
        lambda x: x == x.max())
    logo_df = iteration_df.groupby(['site', 'amino_acid'], as_index=False).agg({
        'relative_likelihood': 'sum'  # Sum relative likelihoods for duplicates
    })
    pwm = logo_df.pivot(index='site', columns='amino_acid', values='relative_likelihood').fillna(0)
    plt.figure(figsize=(max(10, pwm.shape[0]), 4))
    logo = logomaker.Logo(pwm, color_scheme='chemistry')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.ax.set_ylabel('Relative Likelihood')
    logo.ax.set_xlabel('Codon Site')
    plt.title('Amino Acid Logo Plot Based on Relative Likelihoods')
    if chain == "IGH":
            plt.savefig("amino_acid_logo_plot.png", dpi=300, bbox_inches='tight')
            plt.close('all')
    else: 
        plt.savefig("amino_acid_logo_plot_light.png", dpi=300, bbox_inches='tight')
        plt.close('all')
    logo_df = iteration_df.groupby(['nt_site', 'nt'], as_index=False).agg({
        'relative_likelihood': 'sum'})
    pwm_nuc = logo_df.pivot(index='nt_site', columns='nt', values='relative_likelihood').fillna(0)
    plt.figure(figsize=(max(10, pwm_nuc.shape[0]), 4))
    logo = logomaker.Logo(pwm_nuc, color_scheme='chemistry')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.ax.set_ylabel('Relative Likelihood')
    logo.ax.set_xlabel('Site')
    plt.title('Nucleotide Logo Plot Based on Relative Likelihoods')
    if chain == "IGH":
            plt.savefig("nucleotide_logo_plot.png", dpi=300, bbox_inches='tight')
            plt.close('all')
    else: 
        plt.savefig("nucleotide_logo_plot_light.png", dpi=300, bbox_inches='tight')
        plt.close('all')
        
    try:
        if chain == "IGH":
            with open(base_string + "/UCA.txt", "w") as f:
                f.write(new_germline + "\n")
            with open(base_string + "/UCA_lhoods.txt", "w") as f:
                f.write(str(new_lhoods) + "\n")
            data.to_csv(base_string + "/UCA_data.csv", index=False)
            iteration_df.to_csv(base_string + "/recombination_stats.csv", index=False)
            pwm.to_csv(base_string + "/logo_plot_pwm.csv", index=False) 
            pwm_nuc.to_csv(base_string + "/logo_plot_pwm_nuc.csv", index=False)
            collapsed_df.to_csv(base_string + "/collapsed_iteration_df_amino_acid.csv", index=False)
        else:
            with open(base_string + "/UCA_light.txt", "w") as f:
                f.write(new_germline + "\n")
            with open(base_string + "/UCA_lhoods_light.txt", "w") as f:
                f.write(str(new_lhoods) + "\n")
            data.to_csv(base_string + "/UCA_data_light.csv", index=False)
            iteration_df.to_csv(base_string + "/recombination_stats_light.csv", index=False)
            pwm.to_csv(base_string + "/logo_plot_pwm_light.csv", index=False) 
            pwm_nuc.to_csv(base_string + "/logo_plot_pwm_nuc_light.csv", index=False)
            collapsed_df.to_csv(base_string + "/collapsed_iteration_df_amino_acid_light.csv", index=False)

    except Exception as e:
        print(f"Error writing output files: {e}")


if __name__ == '__main__':
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description='Process UCA arguments.')
    parser.add_argument('--args_json', required=True, help='Path to JSON file with all arguments')
    args = parser.parse_args()

    with open(args.args_json) as f:
        params = json.load(f)

    if int(params.get("quiet", 0)) > 0:
        print("Arguments: ", params)

    clone_ids_list = params.get("clone_ids", "").split(',')
    starting_germlines_list = params.get("starting_germlines", "").split(',')
    junction_locations_list = params.get("junction_locations", "").split(',')
    tree_table_list = params.get("tree_tables", "").split(',')
    chains_list = params.get("chains", "").split(',')

    if params.get("search", "") == "codon":
        index_table = pd.DataFrame({
            'clone_ids': clone_ids_list,
            'starting_germline': starting_germlines_list,
            'junction_locations': junction_locations_list,
            'tree_table': tree_table_list,
            'chains': chains_list
        })
        for index, row in index_table.iterrows():
            clone_number = row['clone_ids']
            base_string = params.get("directory", "") + "/" + params.get("id", "") + "_" + clone_number
            if int(params.get("quiet", 0)) > 0:
                print("\nRunning on clone", clone_number, "with locus", row['chains'], "saved at", base_string)
            
            if not os.path.exists(base_string):
                os.makedirs(base_string)
            chain = row['chains']

            table_path = row['tree_table']    
            igphyml_df = pd.read_csv(table_path, sep = "\t", header = None)
            igphyml_df.columns = ["site", "codon", "partial_likelihood", "nope", "nada", "no", "equilibrium"]
            igphyml_df['value'] = igphyml_df['partial_likelihood'] + np.log(igphyml_df['equilibrium'])
            germline_string = row['starting_germline']
            junction_string = row['junction_locations']
            starting_point, ending_point = read_starting_and_ending_points(junction_string)
            with open(germline_string, "r") as f:
                starting_germline = f.read().strip()
            starting_cdr3 = starting_germline[starting_point:ending_point]

            v_gene = starting_germline[0:starting_point]
            j_gene = starting_germline[ending_point:len(starting_germline)]
            if "N" in v_gene or "N" in j_gene:
                if int(params.get("quiet", 0)) > 0:
                    print("N in the J gene -- adding to igphyml df")
                codons = [starting_germline[i:i+3] for i in range(0, len(starting_germline), 3)]
                indices_with_N = [index for index, value in enumerate(codons) if 'N' in value and not (len(v_gene) // 3 + 1 <= index <= len(codons) - len(j_gene) // 3  - 1)]
                new_rows = []
                for index in indices_with_N:
                    sub = igphyml_df[igphyml_df['site'] == index]
                    sub_sum = sub['value'].sum()
                    new_rows.append({'site': index, 'codon': codons[index], 'partial_likelihood': sub_sum, 'nope': 0, 'nada': 0, 'no': 0, 'equilibrium': 0, 'value': 0})
                igphyml_df = pd.concat([igphyml_df, pd.DataFrame(new_rows)], ignore_index=True)
                igphyml_df = igphyml_df.sort_values(by='site', ascending=True)
                # replace the Ns in the cdr3 with C and stitch it back together
                cdr3 = starting_germline.split(v_gene)
                cdr3  = cdr3[1].split(j_gene)[0]
                starting_germline = v_gene + cdr3.replace("N", "C") + j_gene
            
            #igphyml_df['value'] = igphyml_df['partial_likelihood'] + np.log(igphyml_df['equilibrium'])
            
            if "N" in starting_cdr3:
                if int(params.get("quiet", 0)) > 0:
                    print("N in starting CDR3 -- replacing with C")
                starting_cdr3 = starting_cdr3.replace("N", "C")
                # replace the Ns in the cdr3 with C and stitch it back together
                starting_germline = v_gene + starting_cdr3 + j_gene

            # split the germline into a list of codons
            codons = [starting_germline[i:i+3] for i in range(0, len(starting_germline), 3)]

            if chain == "IGH":
                model_folder = params.get("model_folder", "")
            elif chain == "IGK":
                model_folder = params.get("model_folder_igk", "")
            elif chain == "IGL":
                model_folder = params.get("model_folder_igl", "")
            else:
                print("Invalid chain type. Please use IGH, IGK, or IGL.")
                sys.exit(1)

            params_file_name = model_folder + '/model_params.txt'
            marginals_file_name = model_folder + '/model_marginals.txt'
            V_anchor_pos_file = model_folder + '/V_gene_CDR3_anchors.csv'
            J_anchor_pos_file = model_folder + '/J_gene_CDR3_anchors.csv'

            if chain == "IGH":
                genomic_data = load_model.GenomicDataVDJ()
                genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
                generative_model = load_model.GenerativeModelVDJ()
                generative_model.load_and_process_igor_model(marginals_file_name)
                pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
            else:
                genomic_data = load_model.GenomicDataVJ()
                genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
                generative_model = load_model.GenerativeModelVJ()
                generative_model.load_and_process_igor_model(marginals_file_name)
                pgen_model = pgen.GenerationProbabilityVJ(generative_model, genomic_data)

            values = get_updated_germline(starting_germline, starting_cdr3, igphyml_df, pgen_model, starting_point, ending_point, cdr3_only=True, max_iter=int(params.get("max_iters", 100)), nproc=int(params.get("nproc", 1)))

            new_germline = values[0]
            new_lhoods = values[1]
            data = values[2]
            iteration_df = values[3]

            # make the logo plot
            iteration_df['amino_acid'] = iteration_df['codon'].apply(translate_to_amino_acid)
            logo_df = iteration_df.groupby(['site', 'amino_acid'], as_index=False).agg({
                'relative_likelihood': 'sum'  # Sum relative likelihoods for duplicates
            })
            pwm = logo_df.pivot(index='site', columns='amino_acid', values='relative_likelihood').fillna(0)

            plt.figure(figsize=(10, 4))
            logo = logomaker.Logo(pwm, color_scheme='chemistry')
            logo.style_spines(visible=False)
            logo.style_spines(spines=['left', 'bottom'], visible=True)
            logo.ax.set_ylabel('Relative Likelihood')
            logo.ax.set_xlabel('Site')
            plt.title('Amino Acid Logo Plot Based on Relative Likelihoods')

            try:
                if chain == "IGH":
                    with open(base_string + "/UCA.txt", "w") as f:
                        f.write(new_germline + "\n")
                    with open(base_string + "/UCA_lhoods.txt", "w") as f:
                        f.write(str(new_lhoods) + "\n")
                    data.to_csv(base_string + "/UCA_data.csv", index=False)
                    iteration_df.to_csv(base_string + "/recombination_stats.csv", index=False)
                    pwm.to_csv(base_string + "/logo_plot_pwm.csv", index=False) 
                    plt.savefig("amino_acid_logo_plot.png", dpi=300, bbox_inches='tight')
                    plt.close('all')
                else:
                    with open(base_string + "/UCA_light.txt", "w") as f:
                        f.write(new_germline + "\n")
                    with open(base_string + "/UCA_lhoods_light.txt", "w") as f:
                        f.write(str(new_lhoods) + "\n")
                    data.to_csv(base_string + "/UCA_data_light.csv", index=False)
                    iteration_df.to_csv(base_string + "/recombination_stats_light.csv", index=False)
                    pwm.to_csv(base_string + "/logo_plot_pwm_light.csv", index=False) 
                    plt.savefig("amino_acid_logo_plot_light.png", dpi=300, bbox_inches='tight')
                    plt.close('all')
            except Exception as e:
                print(f"Error writing output files: {e}")
    elif params.get("search", "") == "nt":
        index_table = pd.DataFrame({
            'clone_ids': clone_ids_list,
            'starting_germline': starting_germlines_list,
            'junction_locations': junction_locations_list,
            'tree_table': tree_table_list,
            'chains': chains_list,
            'directory': params.get("directory", ""),
            'id': params.get("id", ""),
            'max_iters': params.get("max_iters", 100),
            'quiet': int(params.get("quiet", 0)),
            'model_folder': params.get("model_folder", ""),
            'model_folder_igk': params.get("model_folder_igk", ""),
            'model_folder_igl': params.get("model_folder_igl", "")
        })

        rows = [row.to_dict() for _, row in index_table.iterrows()]
        with Pool(processes=int(params.get("nproc", 1))) as pool:
            pool.map(process_row, rows)
