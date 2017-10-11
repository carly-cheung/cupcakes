# Notes on directories and files

## Commands

### `pwd`

* "print working directory" - shows which folder you're currently in
* like GPS for your computer - tells you where you are

### `cd`

* Changes directory
* Changes the folder you are in
* `cd -` takes you to the *previous* folder, but not the *parent* folder

For example:

```
$ pwd
/Users/olga/Desktop
$ cd ~
$ pwd
/Users/olga
$ cd -
$ pwd
/Users/olga/Desktop
```


### `ls`

- `ls -F`: Lists all files and folders in current directory and changes the formatting by adding a slash.
- `ls -F -a`: In addition to `ls -F`, shows all hidden files
	- Can think of `-a` as "all"
- `ls Desktop`: Looks for a folder in your current directory called `Desktop` and shows you the files that are in `Desktop`. If there's no folder called `Desktop`, it complains.


### `man`

- `man` shows you what the command means.
	- `man` is short for "manual"
- One can exit out of `man` by pressing the letter "q"

## Absolute vs relative paths quiz

1. `cd .`: takes you to the `.` directory - aka where you currently are
2. `cd /`: takes you to the "root" directory - the parent of all folders in your computer
3. `cd /home/amanda`: Doesn't work because there's no folder called `/home` and the home directory starts with the word `Users`