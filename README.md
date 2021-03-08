# Beta-Robust-Scheduling
Codes for Paper "Scheduling Stochastic Jobs - Complexity and Approximation Algorithms" [ICAPS' 2021]

How to run this code?

1. put the execuable file along with the test instance
    
   for example, 
   
   -------------------------
   
   -dataseta
   
   ------0.txt
   
   ------****
   
   ------99.txt
   
   -execuablefile

2. use commond  "execuablefile dataset startno endno config_ip"

   [for example,   execuablefile dataseta 0 100 config_ip]
   
   the optimal answers for instance dataseta/0.txt .... dataseta/99.txt store in the file dataseta/opt_std_config.txt

3. use commond  "execuablefile dataset startno endno approx_fpt valxi numberofthread"

   [for example, execuablefile dataseta 0 100 approx_fpt 10 8]

   the answers of Algorithm 3 (with parameter xi=valxi) for instance dataseta/0.txt .... dataseta/99.txt store in the file dataseta/approx_fpt_valxi.txt

