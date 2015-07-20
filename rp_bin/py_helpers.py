####################
#
# py_helpers.py
# By Raymond Walters, Jul 2015
#
# Set of utility functions for use in picopili
#
####################


# wc -l, taken from http://stackoverflow.com/questions/845058
def file_len(fname):
    import subprocess    
    
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

    

# read ricopili config file as dict    
def read_conf(fname):
    
    configs = {}
    
    with open(fname, 'r') as f:
        for line in f:
            (key, val) = line.split()
            configs[str(key)] = val  
    
    return configs
    


# find mail program
def find_mail():
    # init
    scr_path = None
    from distutils import spawn
    
    # try mutt first
    scr_path = spawn.find_executable("mutt")

    # if not found try mail
    if scr_path == None:
        scr_path = spawn.find_executable("mail")
    
    return scr_path



# send email with given subject using text in file
def send_mail(subj, fname):

        import os
        import subprocess
        
        # get email script
        email_script = find_mail()
        if email_script == None:
            raise IOError("Unable to find 'mutt' or 'mail' in path to send email")
        
        # get mail address from config file
        configs = read_conf(os.environ['HOME']+"/ricopili.conf")

        # verify file before send
        if not os.path.isfile(fname):
            raise IOError("Filename does not exist (%s)" % str(fname))
        
        # send
        subprocess.check_call([email_script,
                               "-s", str(subj),
                               configs['email'],
                               "<", str(fname)])



