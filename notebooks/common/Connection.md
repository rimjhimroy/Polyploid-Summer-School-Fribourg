# Polyploid Summer School 2023
10-20 July 2023

## Connecting to IBU Server with Student IDs and Passwords

Protocol: Connecting to IBU Server with Student IDs and Passwords

1. Requirements:
   - Student ID: The unique identifier assigned to each student.
   - Password: The password associated with the student ID.
   - Server URL: `binfservms01.unibe.ch`

2. Open a Terminal or Command Prompt:
   - On Windows: Press Win+R, type `cmd`, and press Enter.
   - On macOS: Press Cmd+Space, type `Terminal`, and press Enter.

3. Connect to the IBU Server:
   - In the terminal, type the following command:
     ```
     ssh studentXX@binfservms01.unibe.ch
     ```
     Replace `studentXX` with your actual student ID.

4. Enter the Password

5. Successful Connection:
   - If the student ID and password are correct, and the server is accessible, you will now see a command prompt indicating a successful connection.

6. Executing Commands:
   - Once connected to the IBU server, you can execute various commands to interact with the server, such as navigating directories, running programs, or accessing files.
   - Use common Unix/Linux commands like `ls`, `cd`, `mkdir`, `cp`, `mv`, `rm`, etc., to perform operations on the server.

7. Closing the Connection:
   - To close the connection and disconnect from the IBU server, type `exit` and press Enter.
   - You will be returned to your local machine's command prompt.

## Data access and Resource allocation
Please note the following information regarding data access and resource allocation on the cluster:

1. Data Access:
   - The data for the Summer school practical sessions are located in the directory: /data/courses/pploidy/
   - Students are in the Unix group `pploidystudents` have read-only permissions for these data.

2. User Quotas:
   - Each student and user has a 4TB quota to write into the directory /data/users/loginname.
   - Each student and user has a 20GB quota to write into the directory /home/loginname.
   - Use the `lsquota` command to view the individual quota usage.

3. Cluster Resources:
   - The Slurm partition for the cluster is named `pploidy`.
   - The `pploidy` partition has 128 CPUs reserved.
   - The `pploidy` partition has 512 GB of RAM.


## How to submit jobs on the IBU cluster
You can create a job script to submit jobs on the cluster. The job script contains the commands that you would like to run on the cluster. You can submit the job script using the `sbatch` command. The job will be added to the queue and will be executed when the resources are available. You can check the status of the job using the `squeue` or `sacct` commands. You can cancel the job using the `scancel` command.

1. Example on how to create a job submission script (e.g. `submit.sh`) with the following content:
   ```
   #!/bin/bash
   #SBATCH --job-name=polyploid_genotyping
   #SBATCH --output=polyploid_genotyping_%j.out
   #SBATCH --error=polyploid_genotyping_%j.err
   #SBATCH --time=00:30:00
   #SBATCH --mem=10G
   #SBATCH --cpus-per-task=1
   #SBATCH --partition=pploidy

   # Load the miniconda module
   module load Conda/miniconda/latest

   # Activate the Conda environment created for the summer school
   conda activate /data/courses/pploidy/polyploid-genotyping/env/polyploid_env

   # Run the command
   gatk MarkDuplicates -I input.bam -O deduplicated.bam -M metrics.txt
   ```
2. Submit the job using the `sbatch` command:
   ```
   sbatch submit.sh
   ```
3. Check the status of the job using the `squeue` command:
   ```
   squeue -u <your-username>
   ```
4. Check the status of the job using the `sacct` command:

   ```
   sacct -j <job-id>
   ```
5. Cancel the job using the `scancel` command:
   ```
   scancel <job-id>
   ```
Alternatively, you can run the command directly uing the --wrap option:
```
sbatch -J polyploid_genotyping -o polyploid_genotyping_%j.out -e polyploid_genotyping_%j.err --time=00:30:00 --mem=10G -c 1 --wrap="gatk MarkDuplicates -I input.bam -O deduplicated.bam -M metrics.txt"
```