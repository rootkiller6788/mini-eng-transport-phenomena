import os;b=os.path.dirname(os.path.abspath(__file__));
def a(p,c):
    with open(os.path.join(b,p),"a") as fh: fh.write(c)
    return len(c.splitlines())
t=0
