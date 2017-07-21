#coding:utf-8

class ReedSolomonError(Exception):
    pass

def gf_sub(x, y):
    return x ^ y

def gf_mult_noLUT(x, y, prim=0x11d):
    def cl_mult(x, y):
        z = 0
        i = 0
        while (y >> i) > 0:
            if y & (1 << i):
                z ^= x << i
            i += 1
        return z

    def bit_length(n):
        bits = 0
        while n >> bits: bits += 1
        return bits

    def cl_div(dividend, divisor=None):
        dl1 = bit_length(dividend)
        dl2 = bit_length(divisor)
        if dl1 < dl2:
            return dividend
        for i in range(dl1 - dl2, -1, -1):
            if dividend & (1 << i + dl2 - 1):
                dividend ^= divisor << i
        return dividend
    result = cl_mult(x, y)
    if prim > 0:
        result = cl_div(result, prim)

    return result


gf_exp = [0] * 512
gf_log = [0] * 256

def init_tables(prim=0x11d):
    global gf_exp, gf_log
    gf_exp = [0] * 512
    gf_log = [0] * 256
    x = 1
    for i in range(0, 255):
        gf_exp[i] = x
        gf_log[x] = i
        x = gf_mult_noLUT(x, 2, prim)

    for i in range(255, 512):
        gf_exp[i] = gf_exp[i - 255]
    return [gf_log, gf_exp]

def gf_mul(x,y):
    if x==0 or y==0:
        return 0
    return gf_exp[gf_log[x] + gf_log[y]]

def gf_div(x,y):
    if y==0:
        raise ZeroDivisionError()
    if x==0:
        return 0
    return gf_exp[(gf_log[x] + 255 - gf_log[y]) % 255]

def gf_pow(x, power):
    return gf_exp[(gf_log[x] * power) % 255]


def gf_inverse(x):
    return gf_exp[255 - gf_log[x]] # gf_inverse(x) == gf_div(1, x)

def gf_poly_scale(p,x):
    r = [0] * len(p)
    for i in range(0, len(p)):
        r[i] = gf_mul(p[i], x)
    return r

def gf_poly_add(p,q):
    r = [0] * max(len(p),len(q))
    for i in range(0,len(p)):
        r[i+len(r)-len(p)] = p[i]
    for i in range(0,len(q)):
        r[i+len(r)-len(q)] ^= q[i]
    return r

def gf_poly_mul(p,q):
    r = [0] * (len(p)+len(q)-1)
    for j in range(0, len(q)):
        for i in range(0, len(p)):
            r[i+j] ^= gf_mul(p[i], q[j])
    return r

def gf_poly_eval(poly, x):
    if poly!=[]:
        y = poly[0]
        for i in range(1, len(poly)):
            y = gf_mul(y, x) ^ poly[i]
        return y
    return


def rs_generator_poly(nsym):
    g = [1]
    for i in range(0, nsym):
        g = gf_poly_mul(g, [1, gf_pow(2, i)])
        # g=[i%2 for i in g]
    return g


def gf_poly_div(dividend, divisor):
    msg_out = list(dividend)
    for i in range(0, len(dividend) - (len(divisor)-1)):
        coef = msg_out[i]
        if coef != 0:
            for j in range(1, len(divisor)):
                if divisor[j] != 0:
                    msg_out[i + j] ^= gf_mul(divisor[j], coef)
    separator = -(len(divisor)-1)
    return msg_out[:separator], msg_out[separator:]
