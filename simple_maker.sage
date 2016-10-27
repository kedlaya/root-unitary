\load("table_functions.sage")
load("prescribed_roots.sage")
from sage.databases.cremona import cremona_letter_code
from sage.databases.cremona import class_to_int
#from lmfdb.WebNumberField import WebNumberField

#
#For some reason you need to instantiate polyRing.<x> = ZZ[]
#

def make_simples(g,q):
    print "starting g=%s, q=%s..." % (g,q)
    p,r = q.is_prime_power(get_data=True)
    
    #open a file for writing 
    current_file_name = "weil" + "-" + ( "%s" % g ) + "-"+ ("%s" % q) + ".txt"
    target = open(current_file_name, 'w')
    
    polyRing.<x> = PolynomialRing(ZZ)
    F = Qp(p) #p-adic version
    polyRingF.<t> = PolynomialRing(F)
    
    Lpolys,_ = roots_on_unit_circle( 1+ (q*x^2)^g )
    npolys = len(Lpolys)
    for polycount,Lpoly in enumerate(Lpolys):
        if polycount % 500 == 0 and polycount != 0:
            print "%s / %s polynomials complete" % (polycount, npolys)
        #some stuff that will be used later
        Ppoly = Lpoly.reverse()
        coeffs = Lpoly.coefficients(sparse=False)
        Ppolyt = polyRingF(Ppoly) #p-adic (automatically changes variable from x to t)
        
        if Lpoly.is_irreducible():
            invs, newton_slopes = find_invs_and_slopes(p,r,Ppoly)
            e = lcm([a.denominator() for a in invs])
            #check that invs 
            if e==1:
                #return a single factor
                line = '['
                
                #label
                my_label = make_label(g,q,Lpoly)
                line = line + quote_me(my_label) + ','
                
                #dim
                line = line + str(g) + ','
                
                #q
                line = line + str(q) + ','
                
                #polynomial
                line = line + str(coeffs) + ','
                
                #angle numbers
                line = line + str(sorted(angles(Lpoly))) + ','
                
                #angle ranks
                line = line + '"",'
                
                #number_field label
                #line = line + WebNumberField.from_polynomial(Ppoly).label + ','
                
                slopes,p_rank = newton_and_prank(p,r,Ppolyt)
                
                #p-rank
                line = line + str(p_rank) + ','
                
                #slopes
                printed_slopes = '['
                for s in sorted(slopes):
                    if printed_slopes != '[':
                        printed_slopes += ','
                    printed_slopes += quote_me(s)
                #s = [quote_me(x) for x in sorted(slopes)]
                line = line + printed_slopes + '],'
                
                #A_counts
                a_counts = abelian_counts(g,p,r,Lpoly)
                line = line + str(a_counts) + ',' 
                
                #C_counts
                c_counts = curve_counts(g,q,Lpoly)
                line = line + str(c_counts) + ','
                
                #known jacobian
                line = line + '0,'
                
                #principally polarizable
                line = line + '0,'
                
                #decomposition
                line = line + '[' + quote_me(my_label) + ',1],'
                
                #brauer invariants
                printed_invs = '['
                for inv in invs:
                    if printed_invs != '[':
                        printed_invs += ','
                    printed_invs += quote_me(inv)
                #invs_str = [quote_me(inv) for inv in invs] 
                line = line + printed_invs + '],'
                
                #primitive models
                line = line + '""'
                
                target.write(line + ']' + '\n')
            
            #if e is not one we throw this dude away
        else:
            factors = Lpoly.factor()
            if len(factors)==1:
                factor, power = factors[0]
                invs, newton_slopes = find_invs_and_slopes(p,r,factor.reverse())
                e = lcm([a.denominator() for a in invs])
                
                if e == power:
                    #return a single factor
                    
                    line = '['
                    
                    #label
                    my_label = make_label(g,q,Lpoly)
                    line = line + quote_me(my_label) + ','
                    
                    #dim
                    line = line + str(g) + ','
                
                    #q
                    line = line + str(q) + ','
                     
                    #coeffs
                    line = line + str(coeffs) + ','
                    
                    #angle numbers
                    line = line + str(sorted(angles(Lpoly))) + ','
                    
                    #angle ranks
                    line = line + '" ",'
                    
                    #number_field label
                    #line = line + WebNumberField.from_polynomial(Ppoly).label + ','
                    
                    #FIXME THE FOLLOWING IS STUPID
                    slopes,p_rank = newton_and_prank(p,r,Ppolyt)
                    #line = line + ("%s" % slopes) + ','
                    
                    #p-rank
                    line = line + str(p_rank) + ','
                    
                    #slopes
                    printed_slopes = '['
                    for s in sorted(slopes):
                        if printed_slopes != '[':
                            printed_slopes += ','
                        printed_slopes += quote_me(s)
                    #s = [quote_me(x) for x in sorted(slopes)]
                    line = line + printed_slopes + '],'
                    
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
                    line = line + '[' +  quote_me(my_label) + '],'
                    
                    #brauers
                    invs_str = [quote_me(inv) for inv in invs] 
                    line = line + str(invs_str) + ','
                    
                    #primitive models
                    line = line + '""' 
                    
                    target.write(line + ']' + '\n')
            #throw this guy away
