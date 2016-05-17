#non-simples from simples
load("table_functions.sage")

#fix some q and g
def nonsimples_from_simples(filename,Q,M):
    D = load_previous_polys(flename)

q = 5
g = 3
  
for r in g.partition()
    J=dudes_less_than([len(D[ri,q]) for ri in r]):
    dict = {}
    for j in J:
        factors = [ D[ri,q][j[i]] for i in range(len(r))]
        
        Lpoly = product([factor[3] for factor in factors])
        
        #need to compute a method for computing the invs and newton poly
        
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
        s = [quote_me(x) for x in str(sorted(slopes))]
        line = line + str(s) + ','
                
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
        factor_labels = [quote_me(factor[0]) for factor in factors]
        line = line + '[' + quote_me(factor_labels) + '],'
                
        #brauers
        line = line + '" "' + ','
                
        #primitive models
        line = line + '""'
                
        target.write(line + ']' + '\n')
        


#The new L-function is the product of the L_i where the L_i is in the list associated to r_i,q we need to extract the label. 
