from precalculate import *


a=[1,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,1,0,1,1,0,1,1,0,1,0,1,1,0,1]

def rs_encode_msg(msg_in, nsym):
    if (len(msg_in) + nsym) > 255: raise ValueError("Message is too long (%i when max is 255)" % (len(msg_in)+nsym))
    gen = rs_generator_poly(nsym)
    msg_out = [0] * (len(msg_in) + len(gen)-1)
    msg_out[:len(msg_in)] = msg_in

    for i in range(len(msg_in)):
        coef = msg_out[i]
        if coef != 0:
            for j in range(1, len(gen)):
                msg_out[i+j] ^= gf_mul(gen[j], coef) # equivalent to msg_out[i+j] += gf_mul(gen[j], coef)

    msg_out[:len(msg_in)] = msg_in
    return msg_out

# print rs_encode_msg(a,8),'0'

def rs_calc_syndromes(msg, nsym):
    synd = [0] * nsym
    for i in range(0, nsym):
        synd[i] = gf_poly_eval(msg, gf_pow(2,i))
    # synd = [i % 2 for i in synd]
    return [0] + synd
# print rs_calc_syndromes([1,1,0,1,0,1,0],2)
# print rs_calc_syndromes(rs_encode_msg(a,3),8),'1'


def rs_check(msg, nsym):
    return ( max(rs_calc_syndromes(msg, nsym)) == 0 )