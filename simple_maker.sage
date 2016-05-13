load("table_functions.sage")
load("prescribed_roots.sage")
from sage.databases.cremona import cremona_letter_code
from sage.databases.cremona import class_to_int
#from lmfdb.WebNumberField import WebNumberField

def make_simples(g,q):
    p,r = q.is_prime_power(get_data=True)
    
    #open a file for writing 
    current_file_name = "weil" + "-" + ( "%s" % g ) + "-"+ ("%s" % q) + ".txt"
    target = open(current_file_name, 'w')
    
    polyRing.<x> = PolynomialRing(ZZ)
    F = Qp(p) #p-adic version
    polyRingF.<t> = PolynomialRing(F)
    
    Lpolys,_ = roots_on_unit_circle( 1+ (q*x^2)^g )
    
    for Lpoly in Lpolys:
        #some stuff that will be used later
        Ppoly = Lpoly.reverse()
        coeffs = Lpoly.coefficients(sparse=False)
        Ppolyt = polyRingF(Ppoly) #p-adic (automatically changes variable from x to t)
        
        if Lpoly.is_irreducible():
            invs, newton_slopes = find_invs_and_slopes(p,r,Ppoly)
            e = max(invs)
            #check that invs 
            if e==1:
                #return a single factor
                
                #label
                my_label = make_label(g,q,Lpoly)
                line = quote_me(my_label) + ','
                
                #coeffs
                line = line + str(coeffs) + ','
                
                #angle numbers
                line = line + str(angles(Lpoly)) + ','
                
                #number_field label
                #line = line + WebNumberField.from_polynomial(Ppoly).label + ','
                
                #FIXME THE FOLLOWING IS STUPID
                slopes,p_rank = newton_and_prank(p,r,Ppolyt)
                line = line + ("%s" % slopes) + ','
                
                #p-rank
                line = line + str(p_rank) + ','
                
                #slopes
                line = line + str(slopes) + ','
                
                #a-counts
                a_counts = abelian_counts(g,p,r,Lpoly)
                line = line + str(a_counts) + ',' 
                
                #c-counts
                c_counts = curve_counts(g,q,Lpoly)
                line = line + str(c_counts) + ','
                
                #known jacobian
                line = line + '0,'
                
                #principally polarizable
                line = line + '0,'
                
                #decomposition
                line = line + quote_me(my_label) + ','
                
                #brauers
                line = line + str(invs) + ','
                
                #primitive models
                line = line + '""'
                
                target.write(line + '\n')
            
            #if e is not one we throw this dude away
        else:
            factors = Lpoly.factor()
            if len(factors)==1:
                factor, power = factors[0]
                invs, newton_slopes = find_invs_and_slopes(p,r,factor.reverse())
                e = max(invs)
                
                if e == power:
                    #return a single factor
                    
                    #label
                    my_label = make_label(g,q,Lpoly)
                    line = quote_me(my_label) + ','
                     
                    #coeffs
                    line = line + str(coeffs) + ','
                    
                    #angle numbers
                    line = line + str(angles(Lpoly)) + ','
                    
                    #number_field label
                    #line = line + WebNumberField.from_polynomial(Ppoly).label + ','
                    
                    #FIXME THE FOLLOWING IS STUPID
                    slopes,p_rank = newton_and_prank(p,r,Ppolyt)
                    line = line + ("%s" % slopes) + ','
                    
                    #p-rank
                    line = line + str(p_rank) + ','
                    
                    #slopes
                    line = line + str(slopes) + ','
                    
                    #a-counts
                    a_counts = abelian_counts(g,p,r,Lpoly)
                    line = line + str(a_counts) + ',' 
                    
                    #c-counts
                    c_counts = curve_counts(g,q,Lpoly)
                    line = line + str(c_counts) + ','
                    
                    #known jacobian
                    line = line + '0,'
                    
                    #principally polarizable
                    line = line + '0,'
                    
                    #decomposition
                    line = line + quote_me(my_label) + ','
                    
                    #brauers
                    line = line + str(invs) + ','
                    
                    #primitive models
                    line = line + '""' 
                    
                    target.write(line + '\n')
            #throw this guy away
