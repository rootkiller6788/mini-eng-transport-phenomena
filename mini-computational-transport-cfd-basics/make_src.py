import os
base = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(base, 'src')
os.makedirs(src, exist_ok=True)

def write_file(name, lines):
    p = os.path.join(src, name)
    with open(p, 'w') as fh:
        fh.write('
'.join(lines) + '
')
    return len(lines)

