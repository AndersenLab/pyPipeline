#! /user/bin/env python
import subprocess

def md5(filename):
    """Performs an md5 digest on a given file"""
    try:
        subprocess.check_output("which md5sum", shell=True) # Checks to see if md5sum exists.
        md5_command = "md5sum"
    except:
        md5_command = "md5"
    md5_result = subprocess.Popen('%s %s' % (md5_command, filename), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].replace("\n","")
    print md5_result, '%s %s' % (md5_command, filename)
    return md5_result.split("=")[1].strip()
