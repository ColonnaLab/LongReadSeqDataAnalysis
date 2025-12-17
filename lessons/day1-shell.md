### **1. Introducing the Shell**
This page contain a short and adapted version of the Software Carpentry lesson [The shell](https://swcarpentry.github.io/shell-novice/). We will use it as notes to key concepts we will discuss during our lesson. 


##### **Shell and prompt ($)**

- Graphical user interface (GUI) versus command line interface (CLI)
- The most popular Unix shell is Bash (the Bourne Again SHell — so-called because it’s derived from a shell written by [Stephen Bourne](https://en.wikipedia.org/wiki/Stephen_R._Bourne))
- Instructions to computer through shell are given as text in commands (followed by `Enter`) or sequences of commands (scripts) 
- The command line is necessary to interact with remote machines and supercomputers


##### Why genomics uses command-line tools? 

- **Automation** - Handle large numbers of files and repetitive tasks efficiently
- **Computing power** - Access remote servers and high-performance computing clusters
- **Flexibility** - Customize analyses beyond what GUI tools allow (if they even exist)


```diff
+ the prompt symbol $. Prompt might look different  than $. 
user1@vm-corso-colonna:~$
```

Type a command at the prompt (do not type the prompt) and send the command with `Enter`

```diff
user1@vm-corso-colonna:~$ echo focaccia al pomodoro 
focaccia al pomodoro

user1@vm-corso-colonna:~$ ls  
Shell-lesson-data

user1@vm-corso-colonna:~$ ks   
ks: command not found
```


### **2. Navigating Files and Directories**

- File system: the part of the operating system responsible for managing files and directories.
		 It organizes our data into files, which hold information, and directories (also called ‘folders’), which hold files or other directories.

``` diff 
+ COMMAND  pwd stands for print working directory
 ```

user1@vm-corso-colonna:~$ pwd   
/home/user1

**root (/)**

```diff
+ Root (/) is the top-level directory in a Unix/Linux filesystem hierarchy:
+ ~ means home, your space in the filesystem 
+ ~$ means you are at the prompt in your home directory 

user1@vm-corso-colonna:~$ cd /  
user1@vm-corso-colonna:/$ pwd  
user1@vm-corso-colonna:/
│
├── bin/
├── boot/
├── data/
├── dev/
├── etc/
├── home/
+│   └── user1/      # Your home directory (~)
├── lib/
├── lib32/
├── lib64/
├── libx32/
├── lost+found/
├── media/
├── mnt/
├── opt/
├── proc/
├── root/
├── run/
├── sbin/
├── snap/
├── srv/
├── sys/
├── tmp/
├── usr/
└── var/

```

```
user1@vm-corso-colonna:~
│
└─── shell-lesson-data/
     │
     └─── exercise-data/
          ├─── alkanes/         [DIR]
          ├─── animal-counts/   [DIR]
          ├─── creatures/       [DIR]
          ├─── numbers.txt      [FILE]
          ├─── project/         [DIR]
          └─── writing/         [DIR]
```


``` diff 
+ COMMAND ls stands for list directory contents

user1@vm-corso-colonna:~$ ls  
Shell-lesson-data
```


``` diff 
+ COMMAND use --help and man  learning about command options 

user1@vm-corso-colonna:~$  ls --help

user1@vm-corso-colonna:~$  man ls
```

Exploring Other Directories
```
user1@vm-corso-colonna:~$ ls  /home   

/home/
├── ubuntu/        # System user
└── user1-35/      # Course participants
    ├── user1/     # ← You are here
    ├── user2/
    ├── user3/
    ├── ...
    └── user35/
```


``` diff 
+ COMMAND  cd stands for change the shell working directory 

user1@vm-corso-colonna:~$ ls 
shell-lesson-data

user1@vm-corso-colonna:~$ cd shell-lesson-data/

user1@vm-corso-colonna:~/shell-lesson-data$ pwd 
/home/user1/shell-lesson-data
```


**.hidden files** 

>> Go to the [CHALLENGES](https://swcarpentry.github.io/shell-novice/02-filedir.html) of the SC lesson __ls Reading Comprehension__




#####  **2a. General Syntax of a Shell Command**

![syntax](../img/SC-syntax.png)


+ **Options (also referred as switches or flags)** change the behavior of a command  
		single dash (-) short options  
		two dashes (--) long options

+ **Arguments** tell the command what to operate on (e.g. files and directories). 

+ **Parameters:**  options and arguments .

+ *Space* Each part is separated by spaces. If you omit the space between ls and \-F the shell will look for a command called ls-F, which doesn’t exist  
+ *Capitalization*
		ls \-s will display the size of files and directories alongside the names   
		ls \-S will sort the files and directories by size

```diff
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data$ ls -S   
alkanes  animal-counts  creatures  writing  numbers.txt  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data$ ls -s  
total 20  
4 alkanes  4 animal-counts  4 creatures  4 numbers.txt  4 writing


user1@vm-corso-colonna:~/shell-lesson-data$ ls north-pacific-gyre/  
NENE01729A.txt  NENE01736A.txt  NENE01751B.txt  NENE01843A.txt  NENE01971Z.txt  NENE01978B.txt  NENE02040A.txt  NENE02040Z.txt  NENE02043B.txt  goostats.sh  
NENE01729B.txt  NENE01751A.txt  NENE01812A.txt  NENE01843B.txt  NENE01978A.txt  NENE02018B.txt  NENE02040B.txt  NENE02043A.txt  goodiff.sh
```



### **3. Working With Files and Directories**

```diff
! EXERCISE: See where we are with $pwd  and move to  exercise-data/writing, see what it contains 
```

#####  **3a. Create a directory** 

``` diff 
+ COMMAND   mkdir make directories

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$  mkdir thesis

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls   
LittleWomen.txt  haiku.txt  thesis

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls -F thesis/  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$
```

``` diff 
+ COMMAND  the -p option allows mkdir to create a directory with nested subdirectories in a single operation

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ mkdir -p ../project/data ../project/results  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ll
```

``` diff 
+ COMMAND  the -R option to the ls command will list all nested subdirectories within a directory

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls \-FR ../project/  
../project/:

data/  results/  
../project/data:

../project/results:
```

**Good names for files and directories**

1. Don’t use spaces.  
2. Don’t begin the name with \- (dash)  
3. Stick with lowercase letters, numbers, . (period or ‘full stop’), \- (dash) and \_ (underscore)

##### **3b. Create a text file**

```diff
! EXERCISE: See where we are with $pwd  and move to  `thesis` see what it contains 
```

```
+ Text editors are programs used to create, view, and modify plain text files directly from the command line.

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing/thesis$ nano draft.txt
```

* **Control, Ctrl, or ^ Key**  
* **File names and extensions** 


##### **3c. Moving files and directories**

```diff
! EXERCISE: See where we are with $pwv and return to the shell-lesson-data/exercise-data/writing directory  
```

``` diff 
+ COMMAND  mv - move (rename) files

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls thesis/  
draft.txt  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ mv thesis/draft.txt  thesis/quotes.txt  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls thesis/  
quotes.txt
```

``` diff 
+ move quotes.txt into the current working directory “.”

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$  mv thesis/quotes.txt .  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls   
LittleWomen.txt  haiku.txt  quotes.txt  thesis

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls thesis/  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ 
```

##### **3d. Copying files and directories** 

``` diff 
+ COMMAND   cp - copy files and directories

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$  cp quotes.txt thesis/quotations.txt  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls   
LittleWomen.txt  haiku.txt  quotes.txt  thesis  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls thesis/  
quotations.txt  
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ 
```

##### **3e. Removing files and directories**

``` diff 
+ COMMAND rm - remove files or directories

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ rm quotes.txt   
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/writing$ ls  
LittleWomen.txt  haiku.txt  thesis
```

* **Deleting Is Forever**

### **4. Pipes and Filters**

```diff
! EXERCISE: See where we are with $pwv and return to  shell-lesson-data/exercise-data/
```

``` diff 
+ COMMAND wc is the ‘word count’ command: it counts the number of lines, words, and characters in files (returning the values in that order from left to right)

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$  wc cubane.pdb
  20  156 1158 cubane.pdb
```

``` diff 
+ COMMAND the wildcard * is a glob pattern that matches zero or more characters in file and directory names

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$ wc *.pdb 
  20  156 1158 cubane.pdb
  12   84  622 ethane.pdb
   9   57  422 methane.pdb
  30  246 1828 octane.pdb
  21  165 1226 pentane.pdb
  15  111  825 propane.pdb
 107  819 6081 total
```


``` diff 
+ COMMAND the option  -l shows only the number of lines per file

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$ wc -l *.pdb
  20 cubane.pdb
  12 ethane.pdb
   9 methane.pdb
  30 octane.pdb
  21 pentane.pdb
  15 propane.pdb
 107 total
```
##### **4a. Capturing output from commands**

```diff
! EXERCISE: Which of the files in /home/user1/shell-lesson-data/exercise-data/alkanes contains the fewest lines? 
```


``` diff 
+ COMMAND the symbol  > tells the shell to redirect the command’s output to a file instead of printing it to the screen

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$ wc -l *.pdb > lengths.txt
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$ ls 
cubane.pdb  ethane.pdb  lengths.txt  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
```

``` diff 
+ COMMAND  cat gets its name from ‘concatenate’ i.e. join together, and it prints the contents of files one after another. If there is only one file, it just shows what it contains

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$  cat lengths.txt
  20 cubane.pdb
  12 ethane.pdb
   9 methane.pdb
  30 octane.pdb
  21 pentane.pdb
  15 propane.pdb
 107 total
```

``` diff 
+ COMMAND sort is used to  sort lines of text files

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$  sort lengths.txt
   9 methane.pdb
  12 ethane.pdb
  15 propane.pdb
  20 cubane.pdb
  21 pentane.pdb
  30 octane.pdb
 107 total


+ the -r option sorts in reverse order

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$ sort -r  lengths.txt 
 107 total
  30 octane.pdb
  21 pentane.pdb
  20 cubane.pdb
  15 propane.pdb
  12 ethane.pdb
   9 methane.pdb
```

##### **4b Passing output to another command**

``` diff 
+ vertical bar, |, between the two commands is called a pipe. It tells the shell that we want to use the output of the command on the left as the input to the command on the right.

+ COMMAND  head - output the first part of files
+ using the option -n it is possible to specify how many lines we want to print 

user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$ sort  lengths.txt | head -n 3
   9 methane.pdb
  12 ethane.pdb
  15 propane.pdb
```
<img src="../img/piperedirect.png" alt="pipe" style="width:80%; height:auto;">

>> Go to the [CHALLENGES](https://swcarpentry.github.io/shell-novice/04-pipefilter.html) of the SC lesson __Pipe Reading Comprehension__


##### **Shell Scripts** 

```diff
! EXERCISE write previus command in a text file using nano and save it with the name myfirstscript.sh

```

```diff 
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$ nano myfirstscript.sh 
+ save the file using Ctrl-O 
+ exit the text editor using Ctrl-X  

+ run the script: 
user1@vm-corso-colonna:~/shell-lesson-data/exercise-data/alkanes$ bash myfirstscript.sh 

```