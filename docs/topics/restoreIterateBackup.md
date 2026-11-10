**restoreIterateBackup** - *Put back the files that getTimeTreesIterate moved to a backup folder*

Description
--------------------

Assumes no run or BEAST processes are running with this `dir` and `id`.
Deletes any new checkpoint files since the last backup. Restores checkpoint and 
BEAST files from the last or specified backup. The backup folder is deleted.


Usage
--------------------
```
restoreIterateBackup(dir, id, backup = NULL)
```

Arguments
-------------------

dir
:   directory of the analysis

id
:   analysis id

backup
:   backup folder. Default: newest backup of `id` in `dir`.




Value
-------------------

Invisibly, the names of the files and folders that were put back.




See also
-------------------

[getTimeTreesIterate](getTimeTreesIterate.md)






