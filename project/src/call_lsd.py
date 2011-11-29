# To change this template, choose Tools | Templates
# and open the template in the editor.

import os

if __name__ == '__main__':

    mat_cmd = 'matlab >&! -r \"cd(\'/media/DATAPART1/oakland/app/dev-ah/matlab/lsd\'); fid=fopen(\'tmp.txt\',\'w\'); fclose(fid); quit;\"'
    print mat_cmd
    os.system(mat_cmd)
    #os.system('matlab >&! \"quit\"')
    print 'BEAT YA'