import os
import shlex
import subprocess
import datetime

spP = subprocess.Popen
spgo = subprocess.getoutput

# EDNIX_STRICT=1 -> any non-zero return code raises RuntimeError (recommended for validation runs)
STRICT = os.environ.get('EDNIX_STRICT', '0') == '1'


def printcolor(msg, code):
    bcolors = [['HEADER', 'OKBLUE', 'OKCYAN', 'OKGREEN', 'WARNING', 'FAIL', 'ENDC', 'BOLD', 'UNDERLINE'],
               ['\033[95m', '\033[94m', '\033[96m', '\033[92m', '\033[93m', '\033[91m', '\033[0m', '\033[1m', '\033[4m']]
    for i, j in enumerate(bcolors[0]):
        if j == code:
            nl = bcolors[1][i] + str(msg) + bcolors[1][6]
            print(nl)
            return nl


def my_split(s):
    splits0 = s.split(' ')
    splits = []
    param = False
    for j, elem in enumerate(splits0):
        if not param:
            if '"' in elem:
                param = True
                start = j
            else:
                splits.append(elem)
        elif param and '"' in elem:
            param = False
            splits.append(' '.join(splits0[start:j + 1]))
    return splits


def _tokenize(cmd):
    """shlex handles quotes and repeated spaces; falls back to the legacy split on unbalanced quotes."""
    try:
        return shlex.split(cmd)
    except ValueError:
        return [c for c in cmd.split(' ') if c != '']


def _log(diary_name, *chunks):
    with open(diary_name, 'a') as diary:
        for c in chunks:
            diary.write(f'\n{c}')
        diary.write('\n')


def _check_returncode(cmd, returncode, err_txt, diary_name):
    if returncode == 0:
        return
    nl = f'ERROR: command exited with code {returncode}:\n{cmd}\n{err_txt}'
    printcolor(nl, 'FAIL')
    _log(diary_name, nl)
    if STRICT:
        raise RuntimeError(nl)


def run(cmd, diary_name):
    ct = datetime.datetime.now()
    _log(diary_name, ct, cmd)
    proc = spP(_tokenize(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if out:
        printcolor(out.decode('utf-8', errors='replace'), 'ENDC')
        _log(diary_name, out.decode('utf-8', errors='replace'))
    if err:
        printcolor(err.decode('utf-8', errors='replace'), 'WARNING')
        _log(diary_name, err.decode('utf-8', errors='replace'))
    _check_returncode(cmd, proc.returncode, err.decode('utf-8', errors='replace'), diary_name)
    return out, err


def wb(cmd, diary_name):
    ct = datetime.datetime.now()
    _log(diary_name, ct, cmd)
    proc = spP(_tokenize(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if out:
        txt = out.decode('utf-8', errors='replace').split('(')
        if not txt[-1] == 'is the colortable correct?)\n':
            printcolor(out.decode('utf-8', errors='replace'), 'ENDC')
            _log(diary_name, out.decode('utf-8', errors='replace'))
    if err:
        txt = err.decode('utf-8', errors='replace').split(' ')
        if not txt[0] == '\nWARNING:':
            printcolor(err.decode('utf-8', errors='replace'), 'WARNING')
            _log(diary_name, err.decode('utf-8', errors='replace'))
    _check_returncode(cmd, proc.returncode, err.decode('utf-8', errors='replace'), diary_name)
    return out, err


def get(cmd, diary_name):
    ct = datetime.datetime.now()
    _log(diary_name, ct, cmd)
    proc = spP(_tokenize(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if out:
        _log(diary_name, 'Done.')
    if err:
        printcolor(err.decode('utf-8', errors='replace'), 'FAIL')
        _log(diary_name, err.decode('utf-8', errors='replace'))
    _check_returncode(cmd, proc.returncode, err.decode('utf-8', errors='replace'), diary_name)
    return out, err


def do(cmd, diary_name):
    ct = datetime.datetime.now()
    _log(diary_name, ct, cmd)
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    proc = res.stdout.decode('utf-8', errors='replace').rstrip('\n')   # same content as getoutput()
    if proc:
        printcolor(proc, 'ENDC')
        _log(diary_name, proc)
    _check_returncode(cmd, res.returncode, proc, diary_name)
    return proc


def msg(message, diary_name, code='ENDC'):
    ct = datetime.datetime.now()
    printcolor(message, code)
    _log(diary_name, message, ct)


def error(message, diary_name):
    ct = datetime.datetime.now()
    nl = printcolor(message, 'FAIL')
    _log(diary_name, message, ct)
    return nl


def wait(message, diary_name):
    nl = 'Running command:' + message
    msg(nl, diary_name, 'OKGREEN')
    result = subprocess.run(message, shell=True)
    if result.returncode == 0:
        nl_final = 'INFO: Command completed successfully.'
        msg(nl_final, diary_name, 'OKGREEN')
    else:
        nl_final = 'WARNING: Command failed with return code: ' + str(result.returncode)
        msg(nl_final, diary_name, 'WARNING')
    return nl_final
