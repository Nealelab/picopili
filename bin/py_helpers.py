####################
#
# py_helpers.py
# By Raymond Walters, Jul 2015
#
# Set of utility functions for use in picopili
#
####################

# prevent output buffering
def unbuffer_stdout():
    import sys
    import os
    sys.stdout.flush()
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    return 0


# wc -l, taken from http://stackoverflow.com/questions/845058
def file_len(fname):
    import subprocess    
    
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

   
# tail of file (see wc -l above)
def file_tail(fname, n=1):
    import subprocess    
    
    p = subprocess.Popen(['tail', '-n', str(n), fname], stdout=subprocess.PIPE, 
                                                        stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return str(result)


# read picopili config file as dict    
def read_conf(fname):
    
    configs = {}
    
    with open(fname, 'r') as f:
        for line in f:
            # strips '#' comments at end of line
            # otherwise allows aribtrary content (spaces, etc)
            (key,val) = line.split('#',1)[0].rstrip().split(None,1)
            configs[str(key)] = val  
    
    return configs
    


# test executable exists and can be run, prints error or success
def test_exec(fname, name):
    
    import os
    assert os.path.isfile(str(fname)), "%s not found at %s" % (str(name),str(fname))
    assert os.access(str(fname), os.X_OK), "%s not executable (%s)" % (str(name),str(fname))
    print "%s found: %s" % (str(name),str(fname))



# find and verify executables on environment path
def find_from_path(fname, name):
    
    from distutils import spawn
    file_ex = spawn.find_executable(str(fname))
    if file_ex == None:
        raise IOError('Unable to find %s (%s) in search path' % (str(name), str(fname)))
    else:
        print "%s found: %s" % (str(name), str(file_ex))
        return file_ex



# find executables from either config or path
def find_exec(prog, key=None):

    if key is not None:    
        import os
        conffile = os.environ['HOME']+'/picopili.conf'
        
        if os.path.isfile(conffile):
            configs = read_conf(conffile)
            
            if str(key) in configs:
                exloc = configs[str(key)]+'/'+str(prog)
                test_exec(exloc,str(prog))
                return(exloc)
                
            else:
                print "Config file %s is missing extry %s for %s. Will search on path." % (str(conffile),str(key),str(prog))
            
        else:
            print "Failed to find config file %s. Will search for %s on path." % (str(conffile), str(prog))
            
        exloc = find_from_path(str(prog),str(prog))
        test_exec(exloc,str(prog))
        return exloc



# symlink fromfile to tofile and verify
def link(fromfile, tofile, name):
    
    import os
    os.symlink(str(fromfile), str(tofile))
    if not os.path.isfile(str(tofile)):
        raise IOError("Failed to link %s (%s)" % (str(name), str(tofile)) )




# gzip file, wait for success before removing original
def gz_confirm(fname, outname, force):
    
    import subprocess
    outfile = open(str(outname), 'w')
    if force:
        subprocess.check_call(['gzip', '-vfc', str(fname)], stdout=outfile)
    else:
        subprocess.check_call(['gzip', '-vc', str(fname)], stdout=outfile)

    outfile.close()
    subprocess.check_call(['rm', str(fname)])    



# find mail program
def pp_find_mail():
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
def pp_send_mail(subj, fname):

        import os
        import subprocess
        
        # get mail address from config file
        configs = read_conf(os.environ['HOME']+"/picopili.conf")        
        addr = configs['email']
        
        # don't send email 
        if addr is None or '@' not in str(addr):
            print "Email turned off based on config file entry (%s)." % str(addr)
            return 0
        
        # get email script
        email_script = pp_find_mail()
        if email_script == None:
            raise IOError("Unable to find 'mutt' or 'mail' in path to send email")


        # verify file before send
        if not os.path.isfile(fname):
            raise IOError("Filename does not exist (%s)" % str(fname))
        
        # send
        ff = open(fname, 'r')
        subprocess.check_call([email_script,
                               "-s", str(subj),
                               configs['email']],
                               stdin = ff)
        ff.close()




def file_check_email(filename,taskname):
    
    import os
    if os.path.isfile(filename):
        fini_message = '\n\n' + \
                       '##################################################################\n' + \
                       '##### CONGRATULATIONS: \n' + \
                       '##### ' + str(taskname) + ' finished successfully\n' + \
                       '##################################################################\n'
        print fini_message
        with open('success_file', 'w') as suc_file:
            suc_file.write(fini_message)
      
        pp_send_mail(str(taskname)+'_completed', 'success_file')        
        return True
        
    else:
        err_message = '\n\n' + \
                       '##################################################################\n' + \
                       '##### Error: \n' + \
                       '##### ' + str(taskname) + ' failed to complete: \n' + \
                       '##### Final output file ' + filename + ' not found.\n' + \
                       '##################################################################\n'  
        print err_message
        with open('error_file', 'w') as err_file:
            err_file.write(err_message)

        pp_send_mail(str(taskname)+'_failed', 'error_file')        
        return False  


# nicer warning print format
def warn_format(msg, *a):
    return str(msg) + '\n'

