import subprocess
import datetime

spP = subprocess.Popen
spgo = subprocess.getoutput

def printcolor(msg,code):
    bcolors=[['HEADER','OKBLUE','OKCYAN','OKGREEN','WARNING','FAIL','ENDC','BOLD','UNDERLINE'],
             ['\033[95m','\033[94m','\033[96m','\033[92m','\033[93m','\033[91m','\033[0m','\033[1m','\033[4m']]
    for i, j in enumerate(bcolors[0]):
        if j == code:
            nl = bcolors[1][i] + msg + bcolors[1][6]
            print(nl)
            return nl

def my_split(s):

    splits0 = s.split(' ')
    splits = []

    param = False

    # join everything between parenthesis
    for j, elem in enumerate(splits0):
        if not param:
            if '"' in elem:
                param = True
                start = j
            else:
                splits.append(elem)
        elif param and '"' in elem:
            param = False
            splits.append(' '.join(splits0[start:j+1]))
    return splits

def run(cmd,diary_name):
    ct = datetime.datetime.now()
    diary = open(diary_name, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{cmd}')
    cmd2 = cmd.split(' ')
    proc = spP(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if out:
        printcolor(out.decode("utf-8"),'ENDC')
        diary.write(f'\n{out.decode("utf-8")}')
    #else :
    #    diary.write(f'\n So far so good... ')
    if err:
        printcolor(err.decode("utf-8"),'WARNING')
        diary.write(f'\n{err.decode("utf-8")}')
    diary.write(f'\n')
    diary.close()

    return out,err

def wb(cmd,diary_name):
    ct = datetime.datetime.now()
    diary = open(diary_name, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{cmd}')
    cmd2 = cmd.split(' ')
    proc = spP(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if out:
        txt = out.decode("utf-8").split('(')
        if not txt[-1] == 'is the colortable correct?)\n':
            printcolor(out.decode("utf-8"),'ENDC')
            diary.write(f'\n{out.decode("utf-8")}')
    #else :
    #    diary.write(f'\n So far so good... ')
    if err:
        txt = err.decode("utf-8").split(' ')
        if not txt[0]=='\nWARNING:':
            printcolor(err.decode("utf-8"),'WARNING')
            diary.write(f'\n{err.decode("utf-8")}')
    diary.write(f'\n')
    diary.close()

    return out,err

def get(cmd,diary_name):
    ct = datetime.datetime.now()
    diary = open(diary_name, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{cmd}')
    cmd2 = cmd.split(' ')
    proc = spP(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if out:
        diary.write(f'\nDone.')
    #else :
    #    diary.write(f'\n So far so good... ')
    if err:
        printcolor(err.decode("utf-8"), 'FAIL')
        diary.write(f'\n{err.decode("utf-8")}')
    diary.write(f'\n')
    diary.close()

    return out,err

def do(cmd,diary_name):
    ct = datetime.datetime.now()
    diary = open(diary_name, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{cmd}')
    proc = spgo(cmd)
    if proc:
        printcolor(proc,'ENDC')
        diary.write(f'\n{proc}')
    diary.write(f'\n')
    diary.close()
    return proc

def msg(message,diary_name,code):
    ct = datetime.datetime.now()
    diary = open(diary_name, "a")
    printcolor(message,code)
    diary.write(f'\n{message}')
    diary.write(f'\n{ct}')
    diary.write(f'\n')
    diary.close()

def error(message,diary_name):
    ct = datetime.datetime.now()
    diary = open(diary_name, "a")
    nl = printcolor(message, 'FAIL')
    diary.write(f'\n{message}')
    diary.write(f'\n{ct}')
    diary.write(f'\n')
    diary.close()
    return nl

def wait(message,diary_name):
    nl = "Running command:" + message
    msg(nl,diary_name,'OKGREEN')
    result = subprocess.run(message, shell=True)
    if result.returncode == 0:
        nl_final = 'INFO: Command completed successfully.'
        msg(nl_final,diary_name,'OKGREEN')
    else:
        nl_final = 'WARNING: Command failed with return code:' + result.returncode
        msg(nl_final,diary_name,'WARNING')
    return nl_final







