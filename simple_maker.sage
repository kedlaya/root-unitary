load("table_functions.sage")
from sage.databases.cremona import cremona_letter_code
from sage.databases.cremona import class_to_int

make_simples(g,q):
    p,r = q.is_prime_power(get_data=True)
    #open a file for writing 
    current_file_name = "weil" + "-" + ( "%s" % g ) + "-"+ ("%s" % q) + ".txt"
    
    target = open(currentfilename, 'w')
    polyRing.<x> = PolynomialRing(ZZ)
    F = Qp(p) #p-adic version
    polyRingF.<t> = PolynomialRing(F)
    
    LPolys,_ = roots_on_unit_circle( 1+ (q*x^2)^g )
    
    for Lpoly in Lpolys:
        #some stuff that will be used later
        Ppoly = Lpoly.reverse()
        coeffs = Lpoly.coefficients(sparse=False)'
        
        if LPoly.irreducible():
            invs, newton_slopes = find_invs_and_slopes():
            e = max(invs)
            #check that invs 
            if e==1:
                
                
                #return a single factor
            
            #if e is not one we throw this dude away
        else:
            factors = LPoly.factor()
            if len(factors)==1:
                factor, power = factors[0]
                invs, newton_slopes = find_invs_and_slopes(factor):
                e = max(invs)
                
                if e == power:
                     write_entry() 
                
                #if e is not the power we throw this dude away
                
            #throw this guy away
            
            
    
    
            
            
           
