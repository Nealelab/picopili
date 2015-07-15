import subprocess

# wc -l, taken from http://stackoverflow.com/questions/845058
def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])
    
    
def read_conf(fname):
    
    configs = {}
    
    with open(fname, 'r') as f:
        for line in f:
            (key, val) = line.split()
            configs[str(key)] = val  
    
    return configs
    


