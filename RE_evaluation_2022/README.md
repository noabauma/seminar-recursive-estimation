# Grading Recurisve Estimation Programming Exercise


1. Copy all submissions to the `Submissions` folder

2. Copy the `RE_evaluation_2022` folder to Euler with

       scp -r RE_evaluation_2022/ username@euler.ethz.ch:/cluster/scratch/username

   where `username` is your ETH username.

3. Log in to Euler and move to the RE directory

       ssh username@euler.ethz.ch
       cd $SCRATCH/RE_evaluation_2022/

4. Load the MATLAB module with

       module load new matlab/R2019b

5. Submit the EKF grading job with

       bsub -n 48 -R "rusage[mem=2560]" matlab -nodisplay -singleCompThread -r evaluate_ekf

6. Submit the PF grading job with

       bsub -n 48 -R "rusage[mem=2560]" matlab -nodisplay -singleCompThread -r evaluate_pf

7. Once the job is finished copy the logs and results back to your local machine

       cd RE_evaluation_2022/
       scp username@euler.ethz.ch:/cluster/scratch/username/RE_evaluation_2022/\{*.txt,*.mat,lsf*,*.xls\} .