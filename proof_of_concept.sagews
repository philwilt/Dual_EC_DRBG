"""

Phillip M Wilt
https://github.com/phillwilt/Dual_EC_DRBG

Proof of concept based on an attack from Dan Shumow (Microsoft)

Dan's Presentation: http://rump2007.cr.yp.to/15-shumow.pdf

Insert into a worksheet at: http://cloud.sagemath.org

"""
from multiprocessing import Process;

#bitmask
mask = (2**(30*8) - 1)

def gen_candidates(output):
    """
    Given a 30-byte block of output generate a list of potential 32-byte numbers
    """
    print "Generating Candidates";
    S = [];
    max_i = 2**16; 
    
    #brute force find missing 2 bytes
    for i in range(0,max_i): 
        sh = i << (30*8);         
        x = sh | output;
        z = Mod(x^3-3*x + b, prime256);
        if(z.is_square()):
            y = z.sqrt();
            S.append([x,y]);
    return S;


def is_on_curve(point, curve):
    """
    Checks to see if the point is on the curve. 
    """
    try:
        p = curve(point);
        return true;
    except:
        return false;
        

def test_match(e,Q,curve,pt_list,next_output):
    """
    Given next 2 bytes of output from a Dual_EC_PRNG instance, find the internal state
    
    Params:
        e - d^1 (mod p_ord)
        Q - e*Q = P
        curve - the curve
        pt_list - list of potential points
        output - next 2 bytes of output from a Dual_EC_PRNG instance
    """
    print "Testing matches: " + str(len(pt_list)) + '\n';

    for p in pt_list:
        if(is_on_curve(p, E)):
            A = curve(p); # A = t = r*Q

            s1 = (e*A)[0].lift(); # x-coord: e*A = e*(r0*Q) = e*(r0*d*P) = r0*P
            r1 = (s1*P)[0].lift(); # x-coord: s1 * P = r1
            t1 = (r1*Q)[0].lift(); #x-coord: r1 * Q
            
            pred = t1 & mask; #throw away 16 Most Sig Bits
            test_pred = pred >> (28*8); #get 2 Most Sig 
            
            if(test_pred == next_output):
                print "Match: ",pred;
                print "A: ", A
    
def predict_next(output,next_output, curve):
    """
    Given output, finds internal state of a Dual_EC_PRNG
    
    output - 30-byte output block 
    next_output - next 30-byte output block 
    curve - the curve
    """
    # we only need 2 bytes of Dual_EC_PRNG output to determine state
    output_test = (next_output >> (28*8)) 
    
    #Brute Force Generate Missing Point Data
    print "Brute Force Generation";
    pt_list = gen_candidates(output);
    print "Total Potential Matches: ", len(pt_list);
  
    procs = [];
    jump = 1000;
    
    #Test for matches
    for i in range(0, len(pt_list), jump):
        proc = Process(target=test_match, args=(e,Q,curve,pt_list[i:i + jump],output_test))
        procs.append(proc);
        proc.start();
    
        
    for proc in procs:
        if(proc.is_alive()):
            proc.join()
    


def dual_ec_gen(curve, P, Q, length, seed=None):
    """
    Dual_EC_DRNG instance
    
    Params:
        curve - the curve
        P - point P
        Q - point Q
        length - quantity of numbers to generate
        seed - initial seed (could be an internal state of a comprimised instance)
    """
    rand_list = [];
    # random initial seed
    if(seed == None):
        rand = floor((2**16-1)*random());
        s = int(rand); 
    else:
        s = seed;
       
    for i in range(length):
            r = (s*P)[0].lift(); 
            s = (r*P)[0].lift();
            t = (r*Q)[0].lift();

            rand = t & mask;
            rand_list.append(rand);
    
    return rand_list;


#Curve P-256 from NIST SP800-90
prime256 = 115792089210356248762697446949407573530086143415290314195533631308867097853951;
b = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b; 
E = EllipticCurve(GF(prime256), [0,0,0,-3,b]); #NIST y^2= x^3- 3x + b (mod p)
 
Px = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296; 
Py = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5;
P = E(Px,Py);  

#Original Q
#Qx = 0xc97445f45cdef9f0d3e05e1e585fc297235b82b5be8ff3efca67c59852018192;
#Qy = 0xb28ef557ba31dfcbdd21ac46e2a91e3c304f44cb87058ada2cb815151e610046;

#Backdoor Q
d = 13;
Q = d*P;

p_ord = P.additive_order(); 

e = inverse_mod(d, p_ord); # find e = d^-1 (mod p_ord)
print "e*d: ", Mod(e*d,p_ord); # 1 

print "P: ", P
print "e*Q: ", e*Q

# Spin up a generator
rands = dual_ec_gen(E,P,Q,4);
print "Dual_EC_PRNG Output: "
print rands;

# Drop the hammer
%time q = predict_next(rands[1], rands[2], E);









