import math

def approx_cmp(a,b, approx=True):
    if approx:
        if abs(a-b) > abs(a)*pow(10,-6):
            return True
        else 
            return False
    else:
        if a != b:
            return True
        else
            return False
