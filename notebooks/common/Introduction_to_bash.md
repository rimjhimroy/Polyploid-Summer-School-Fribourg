# Bash Tutorial
Polyploid Summer School 2023
10 - 20 July 2023

## Introduction to Bash
Bash s a popular command-line shell and scripting language used in Linux and macOS systems. It is an essential skill for bioinformatics students as it provides powerful tools for data manipulation, file processing, and automation. In this tutorial, we will cover the basics of Bash scripting with a focus on bioinformatics-related tasks.

## Table of Contents
1. [Getting Started](#getting-started)
   - [Opening the Terminal](#opening-the-terminal)
   - [Basic Commands](#basic-commands)
2. [Working with Files and Directories](#working-with-files-and-directories)
   - [Navigating the Filesystem](#navigating-the-filesystem)
   - [Creating and Removing Files/Directories](#creating-and-removing-filesdirectories)
   - [Copying, Moving, and Renaming Files/Directories](#copying-moving-and-renaming-filesdirectories)
3. [Text Processing](#text-processing)
   - [Viewing File Contents](#viewing-file-contents)
   - [Searching and Filtering Data](#searching-and-filtering-data)
   - [Modifying File Contents](#modifying-file-contents)
4. [Bioinformatics Tasks](#bioinformatics-tasks)
   - [Working with FASTA Files](#working-with-fasta-files)
   - [Counting Sequences](#counting-sequences)
   - [Extracting Sequence Information](#extracting-sequence-information)
5. [Shell Scripting](#shell-scripting)
   - [Creating and Running a Script](#creating-and-running-a-script)
   - [Passing Arguments to a Script](#passing-arguments-to-a-script)
6. [Conclusion](#conclusion)

## Getting Started <a name="getting-started"></a>

### Opening the Terminal <a name="opening-the-terminal"></a>
To start using Bash, open the terminal application on your Linux or macOS system. You can usually find it in the "Utilities" folder or search for it using the system's search function.

### Basic Commands <a name="basic-commands"></a>
Here are some essential Bash commands to get you started:

- `ls` - List files and directories in the current directory.
- `cd` - Change directory.
- `pwd` - Print the current working directory.
- `mkdir` - Create a new directory.
- `rm` - Remove a file or directory.
- `cp` - Copy files and directories.
- `mv` - Move or rename files and directories.
- `cat` - Concatenate and display file contents.
- `grep` - Search for patterns in files.
- `sed` - Stream editor for modifying file contents.
- `head` - Display the first few lines of a file.
- `tail` - Display the last few lines of a file.

## Working with Files and Directories <a name="working-with-files-and-directories"></a>

### Navigating the Filesystem <a name="navigating-the-filesystem"></a>
To navigate the filesystem, you can use the `cd` command followed by the directory path. Here are some examples:

- `cd` - Change to the home directory.
- `cd /path/to/directory` - Change to a specific directory.
- `cd ..` - Move up one directory level.
- `cd ~username` - Change to the home directory of a specific user.

### Creating and Removing Files/Directories <a name="creating-and-removing-filesdirectories"></a>
To create files and directories, you can use the `mkdir` and `touch` commands:

- `mkdir directory_name` - Create a new directory.
- `touch file_name` - Create a new empty file.

To remove files and directories, use the `rm` command:

- `rm file_name` - Remove a file.
- `rm -r directory_name` - Remove a directory and its contents recursively.

### Copying, Moving, and Renaming Files/Directories <a name="copying-moving-and-renaming-filesdirectories"></a>
To copy files or directories, use the `cp` command:

- `cp source_file destination_file` - Copy a file.
- `cp -r source_directory destination_directory` - Copy a directory and its contents recursively.

To move or rename files and directories, use the `mv` command:

- `mv old_name new_name` - Rename a file or directory.
- `mv file_name directory_name` - Move a file to a different directory.

## Text Processing <a name="text-processing"></a>

### Viewing File Contents <a name="viewing-file-contents"></a>
To view the contents of a file, you can use the `cat`, `head`, or `tail` commands:

- `cat file_name` - Display the entire file.
- `head -n 10 file_name` - Display the first 10 lines of a file.
- `tail -n 5 file_name` - Display the last 5 lines of a file.

### Searching and Filtering Data <a name="searching-and-filtering-data"></a>
To search for patterns in files, use the `grep` command:

- `grep pattern file_name` - Search for a pattern in a file.
- `grep -i pattern file_name` - Perform a case-insensitive search.
- `grep -r pattern directory_name` - Search for a pattern recursively in a directory.

### Modifying File Contents <a name="modifying-file-contents"></a>
To modify file contents, you can use the `sed` command:

- `sed 's/pattern/replacement/g' file_name` - Substitute a pattern with a replacement in a file.

## Bioinformatics Tasks <a name="bioinformatics-tasks"></a>

### Working with FASTA Files <a name="working-with-fasta-files"></a>
FASTA files are commonly used to store DNA or protein sequences. Here's an example of working with FASTA files using Bash commands:

- Extract sequences from a FASTA file using `grep`:
  ```bash
  grep -v ">" sequences.fasta > clean_sequences.fasta
  ```

### Counting Sequences <a name="counting-sequences"></a>
To count the number of sequences in a FASTA file, you can use Bash commands along with `grep` and `wc`:

```bash
grep -c ">" sequences.fasta
```

### Extracting Sequence Information <a name="extracting-sequence-information"></a>
To extract specific sequence information from a file, you can use tools like `awk` or `cut`. For example, to extract the sequence lengths from a FASTA file:

```bash
grep -v ">" sequences.fasta | awk '{ print length($0) }'
```

## Shell Scripting <a name="shell-scripting"></a>

### Creating and Running a Script <a name="creating-and-running-a-script"></a>
To create a shell script, you can use a text editor to write a series of Bash commands and save the file with a `.sh` extension. For example:

```bash
#!/bin/bash

# This is a simple script
echo "Hello, World!"
```

Save the script as `hello.sh`. To run the script, use the following command:

```bash
bash hello.sh
```

### Passing Arguments to a Script <a

 name="passing-arguments-to-a-script"></a>
You can pass arguments to a shell script by referencing them using the `$` symbol. For example, modify the `hello.sh` script to accept a name argument:

```bash
#!/bin/bash

# This is a script that greets a person
echo "Hello, $1!"
```

Run the modified script and provide a name argument:

```bash
bash hello.sh John
```

The output will be: "Hello, John!"

## Conclusion <a name="conclusion"></a>

In this tutorial, we covered the basics of Bash scripting for bioinformatics students. You learned how to navigate the filesystem, work with files and directories, process text, perform common bioinformatics tasks, and create shell scripts. Bash provides a wide range of tools and possibilities for data manipulation and automation, empowering you in your bioinformatics work.

Explore further Bash resources, practice writing scripts, and experiment with different commands to expand your knowledge and efficiency in bioinformatics tasks.

Happy scripting!

