import sys
import itertools
import random
from math import floor

#Details.

#author = Dylan Johnston
#description = "GL3(F_p) Representations Toolkit (incls littlewood-richardson, decomp into irreps (in groth group), Steinberg's theorem"
#url =(none yet)
#python required -> 3.6

#last updated: 3/9/19 23:00

#########################################################################################################
#### TO USE: Just compile and you will be asked for inputs in the shell, no need to do anything here ####
#########################################################################################################

#List of methods:

#reflectH1 - reflects in H(1_n) for some specified n
#reflectH2
#reflectH3
#W_dim - weyl dimension formula
#det 
#schur_num - numerator of schur poly,used in s()
#s - evals schur poly
#upto - gives three values n_i which represent the max n_i such that the hyperplane H_(i,n_i) lives "to the left/below" the given weight.
#W_to_F 
#F_to_W
#SLP - Strong Linkage Principle, note this doesnt give you the co-effs, just the terms that could/do appear.
#Steinberg - performs steinbergs:
#coeff_checker - used by the coeff finders to check a set of coeffs (bwloe two are done via iterating through lots of possible combos of coeffs
#SLP_coeff_finder - given a list from SLP() it finds the coeffs
#LR_coeff_finder - given WxW = bunch of Ws it finds the coeffs
#LittleRich - performs littlewood-ricardson completely.

#options_menu - gives options of things you can do in the toolkit
#retry - doesnt recognise input, lets you retry
#thankyou - a final thonk you for using the toolkit
#finalask - asks if you want to go again or end.
#options1-5 - these call on the above math-based methods to give answers.
#main 

def reflectH1(aa,n,p):  #reflects weight 'aa' in hyperplane H_(1,n) for prime p
    return [aa[0],aa[2]+n*p-1,aa[1]+1-n*p]

def reflectH2(aa,n,p):  #reflects weight in hyperplane H_(2,n) for prime p
    return [aa[2]+n*p-2,aa[1],aa[0]+2-n*p]

def reflectH3(aa,n,p):  #reflects weight in hyperplane H_(3,n) for prime p
    return [aa[1]+n*p-1,aa[0]+1-n*p,aa[2]]

def W_dim(aa): #Weyl dimension formula for W(a,b,c)
    return 1/2*(aa[0]-aa[1]+1)*(aa[1]-aa[2]+1)*(aa[0]-aa[2]+2)

def det(M =[]):
    return M[0]*(M[4]*M[8]-M[5]*M[7]) - M[1]*(M[3]*M[8]-M[5]*M[6]) + M[2]*(M[3]*M[7]-M[4]*M[6])

def schur_num(a,b,c,x,y,z,s=1): #numerator formula for schur poly (get denom formula by setting a=b=c=0). This comes from Jacobi's bialternant formula.
    S = [x for x in range(9)]
    S[0],S[1],S[2],S[3] = x**(s*(a+2)), y**(s*(a+2)), z**(s*(a+2)),x**(s*(b+1))
    S[4],S[5],S[6],S[7],S[8] = y**(s*(b+1)), z**(s*(b+1)),x**(s*c), y**(s*c), z**(s*c)    
    return det(S)

rand1,rand2,rand3 = round(random.uniform(4,2),2),round(random.uniform(4,2),2),round(random.uniform(4,2),2) #some rand values to evaluate s on by default.
def s(a,b,c,t=1,v1=rand1,v2=rand2,v3=rand3): #this is the evaluation of the schur poly, t stands for twist.
    p = schur_num(a,b,c,v1,v2,v3,t)    
    q = schur_num(0,0,0,v1,v2,v3,t)
    return p/q

def upto(aa,p):  #this gives a max value 'up_to(i)' for each i=1,2,3 such that all the hyperplanes H_(i,n), n <= up_to(i) are left/south of the point aa for a prime p
    up_to1=0
    up_to2=0
    up_to3=0
    while aa[1]-aa[2] >(up_to1+1)*p - 1:
        up_to1 += 1 
    while aa[0]-aa[2]>(up_to2 +1)*p - 2:
        up_to2 += 1
    while aa[0]-aa[1]>(up_to3 +1)*p - 1:
        up_to3 += 1
    return [up_to1,up_to2,up_to3]


def W_to_F(aa,p):  #decomposes W(aa) into sum of irreducible F terms for aa in the p-restricted region ONLY.
    
    ans = [aa]  #no matter what F(aa) appears in decomp of W(aa) so put it in answer now
    
    if 0 <= aa[0] - aa[1] < p and 0 <= aa[0] - aa[1] < p:  #check if we're in p-restricted region
        if aa[0] - aa[2] <= p-2 or aa[0] - aa[1] == p-1 or aa[1] - aa[2] == p-1:  #if this is true we're in the lower alcove 'or' boundary 
            pass  #do nothing, this is just in here for my own notion of completenesss          
        else:  #W(a,b,c) must be in upper alcove.
            ans.append(reflectH2(aa,1,p))  #we must add the reflection to ans, its the reflection in H_(2,1).
       
        return ans  #returns weights of components as a list of two weights. 
              
    else:  #we're not in p-restricted region
        return "error, not in p-restricted region" #this will never happen when W_to_F() is called from other methods, this condition will be checked beforehand.

def F_to_W(ff,p): ##ff is our weight (a,b,c) belonging to F(a,b,c), p is our prime as usual
    
    ans = [ff] ##ff is in the answer no matter what so may as well put it in the list to be returned now.

    if 0 <= ff[0] - ff[1] < p and 0 <= ff[0] - ff[1] < p: #check if we're in p-restricted region        
        if ff[0] - ff[2] <= p-2 or ff[0] - ff[1] == p-1 or ff[1] - ff[2] == p-1: #lower alcove 'or' boundary conditions
            ans = ans #do nothing, this is just in here for my own notion of completenesss          
        else: #(a,b,c) must be in upper alcove.
            ans.append(reflectH2(ff,1,p)) #we must add the reflection to ans, REMEMBER this has a minus in front of it when printing/doing calcs.
            
    return ans #ans returned as a list of weights
    

def SLP(aa,p): #aa refers to the aa in W(aa) which you want to decompose. Note: This works for W in p-res.

    possibleweights = [aa] ##aa always appears so just include it straight away, also keeps it in first position in the list which is nice for later.
    
    for i in range(0,upto(aa,p)[1]): #this gives us a "beam" running from the point towards the origin in a SW direction perpendicular to the H_(2,n)
        r2 =  reflectH2(possibleweights[i],upto(possibleweights[i],p)[1],p) #reflect each weight in nearest hyperplane H_2,i 
        possibleweights.append(r2)

    for xx in possibleweights: 
        while True: #just keep looping, the breaks within the if statements will get us out when we need to.
            #beam out in South direction first, we'll do NW after (same code pretty much)
            i=0
            bbS = [xx] #bb stands for points obtained in the beam-out south, xx is in there trivially (don't reflect).
            r1 = reflectH1(bbS[i],upto(bbS[i],p)[0],p)
            if r1[0]+r1[1] + 2*r1[2] < 0: #when the refection goes too far down so it cant possibly produce more new terms in p-restricted region this cond. is true.
                 break
            equal = 0 #resets it for each new r being tested
            for yy in possibleweights:                
                match1st = yy[0] - r1[0]
                match2nd = yy[1] - r1[1] #dont need to test 3rd since sum of points is constant so if first 2 match the 3rd will also
                if match1st ==0 and match2nd == 0:
                    equal += 1
            if equal == 0: #only if rr isnt already in possible weights will we add it otherwise this beam-out would never end.
                bbS.append(r1)
                possibleweights.append(r1)
                i +=1
            else:
                break #if rr in already in possible weights then theres no need to do anymore reflections as they will get caught
                      #when rr which is already in possweights gets it's turn in the 'for xx in...' loop.
            #NW direction beamout
            j=0 
            bbNW = [xx] #bb stands for points obtained in the NW beam-out
            r3 = reflectH3(bbNW[j],upto(bbNW[j],p)[2],p)
            if r3[1]+r3[2]-2*r3[0] > 0: #cond when the refection goes too leftwards so it cant possibly produce more new terms in p-restricted region.
                 break
            equal = 0 #resets it
            for yy in possibleweights: #compares r to every weight currently in 'possibleweights' to see if it's new.          
                match1st = yy[0] - r3[0]
                match2nd = yy[1] - r3[1] #dont need to test 3rd since sum of points is constant so if first 2 match the 3rd will also
                if match1st ==0 and match2nd == 0:
                    equal += 1
            if equal == 0: #only if r isnt already in possible weights will we add it otherwise this beamout would never end.
                bbNW.append(r3)
                possibleweights.append(r3)
                j +=1
            else:#if r is already in poss weights then there's no need to do anymore reflects for this xx as they will get found by the same term already in poss wghts later
                break
            
    for pp in possibleweights[:]:
        if pp[0]< pp[1] or pp[1] < pp[2]: #condition to remove xx that aren't in the the p-restricted region
            possibleweights.remove(pp)
    return possibleweights #a list of lists.

def Steinberg(aa,p): 

    upto1,upto2,upto3 = upto(aa,p)[0],upto(aa,p)[1],upto(aa,p)[2]
    
    #element (n,n,0)^(some power) \in p-res region needed if upto = [n,n,0] or [n,n+1,0] depending on rhombus R_{i,j} position, each [1,1,0] or [1,2,0] represents shifting p units up
    #element (n,0,0)^(some power) in p res is needed if upto = [0,n,n] or [0,n+1,n] depending on rhombus position, each [0,1,1] or [,2,1] represents shifting p units right
    #we can reach any triangle/rhombus via a linear combo of up and right moves starting in the main bi-alcove
    #so we see that upto1 = #(1,1,0) basis needed and upto3 = #(1,0,0) needed in powered tensors.
    #we also deduce if F(a,b,c) lives in p^2-restricted region we only need 2 tensors
    #if F(a,b,c) lies outside p^2-restricted region but in the p^3 restricted region we need 3 tensors etc etc.

    n = 1 #this is (will be) the 'n' in the statement: "F(a,b,c) lives in the p**n restricted region but NOT the p**(n-1) restricted region".
    while aa[0]-aa[1] >= p**n or aa[1]-aa[2] >= p**n:
        n+=1

    products=[] #to be a list of all the tensors involved in the product
    aval = aa[0]
    bval = aa[1]
    cval = aa[2]
    
    for i in range(n,0,-1): #decreasing list for n down to 0 (inclusive)
        if i == 1: #if we've reached the p-restricted region with our remaining values just add them to the product of tensors and we're done.
            products.append([aval,bval,cval])
            break
        firstbasis=0 #These are the number of (1,1,0) and (1,0,0) terms appearing in each F term  in the tensor decomposition
        secondbasis=0 #They need to be reset every time.
        while upto1 >= floor(p**(i-2)): #Here we are "scaling" back F(a,b,c) to find its relative position in the p-restricted region. This is just translations.
            firstbasis += 1
            upto1 -= p**(i-2)
            
        while upto3 >= floor(p**(i-2)):
            secondbasis += 1
            upto3 -= p**(i-2)
        aval -= p**(i-1)*(firstbasis+secondbasis)
        bval -= p**(i-1)*firstbasis
        products.append([firstbasis+secondbasis,firstbasis,0]) #adding 1stbasis(1,1,0) + 2ndbasis(1,0,0)
                   
    return products #list of lists, there's been no account of the power here, just remember to calculate with it later.

def coeff_checker(total, evals, coeffs): #total is s(a,b,c)(x,y,z) , evals is a list of schur evals at (x,y,z), coeffs are a list of potential coeffs.
    #this will be called by both strong linkage and L-R rule
    e = 0.00001
    hold = total
    for i in range(len(coeffs)):
        hold = hold - evals[i]*coeffs[i]
    if abs(hold) < e:
        return True
    return False
        
def SLP_coeff_finder(LHS,RHS,p):  #LHS is [a,b,c], RHS is list of W(weights) based on steinburgs and expressing Fs as Ws,p prime
                    
    max_value = int(W_dim(LHS)) #a very loose bound, moreso just so we have one to put into next part, we'll never come nere to hitting it.     
    for k in range(max_value + 1): #increases the bound on the range on each entry, this means we dont check silly things like (0,0,...,200) before (1,1,...1)s
        ranges = [] #empty it
        for i in range(len(RHS)):  
            if k <= max_value:
                ranges.append(range(0,k))
            else:
                ranges.append(range(0,max_value))                   
    
        if k==1: #this is just xx = (0,0,..0) which obviously isnt correct (dim =0), so we'll use this chance to check (1,1,...1) early since it comes up so often
                ranges =[] #empty it, it'll refill for k=2 anyway so no worries there
                for i in range(len(RHS)):
                    ranges.append(range(1,2)) #just fill it with ones.
        for xx in itertools.product(*ranges):
            successes = 0 #define a counter which resets for each xx
            for t in range(len(RHS)): #we will check each set of coeffs at 'length(RHS)' sets of evals of the schur functions. Creating a nxn linear system (in general)
                r1 = round(random.uniform(0.3,2),3)
                r2 = 0
                r3 = 0
                while r2 == 0: #if the random numbers are not distinct it all breaks down since we divide by det above in our definition, so need to do this
                    hold2 = round(random.uniform(0.3,2),3)
                    if hold2 != r1:
                        r2 = hold2
                while r3 == 0:
                    hold3 = round(random.uniform(0.3,2),3)
                    if hold3 != r1 and hold3 != r2:
                        r3 = hold3
                LHSvalue = s(LHS[0],LHS[1],LHS[2],1,r1,r2,r3)
                
                ##the next two blocks deal with calculating RHS at the given random numbers
                inner = [] ##this block subtracts the innermost lists together (this is to deal with any W - W terms that may have appeared from expressing Fs as Ws
                for i in range(len(RHS)):
                    holder_list=[] #initalise/reset the holder.
                    for j in range(len(RHS[i])):
                        power = p**(len(RHS[i]) -1 -j) #this determines if schu function gets evaluated at powers of p instead bc of algebraic rep
                        subtraction = 0
                        for k in range(len(RHS[i][j])):
                            subtraction += (-1)**k*s(RHS[i][j][k][0],RHS[i][j][k][1],RHS[i][j][k][2],1,r1**power,r2**power,r3**power) #deals with fact it's algebraic reps
                        holder_list.append(subtraction)
                    inner.append(holder_list)

                RHSvalues = [] ##this deals with the multiplication part of the nested lists, easy as all the work had to done above
                for i in range(len(inner)):
                    multi = 1
                    for j in range(len(inner[i])):
                        multi *= inner[i][j]
                    RHSvalues.append(multi)
                    
                if coeff_checker(LHSvalue,RHSvalues,xx) == True: #if it works for the random set of values, add 1 to successes and go again.
                    successes += 1
                if coeff_checker(LHSvalue,RHSvalues,xx) == False: #if it doesnt work for a set of rand vals we can scrap the coeffs as they dont work, no need to keep checking.
                    break
            if successes == len(RHS): #if coeffs worked for every set of random numbers we tested it at then return those coeffs as the mostly certainly correct ones. 
                return xx #returns a list of coeffs (ints). List has same length as RHS of course.
    #It didnt work :( It must have been a one in a million round-off error so we'll try again"
    return SLP_coeff_finder(LHS,RHS,p)

def LR_coeff_finder(LHS,RHS,bound):  #LHS is [aa,bb], RHS is list of possible terms (list of lists), bound is just the bound on the coeffs found by Weyl-dim
    #this function is super similar to above, unfortunately it differs slightly so I found it easier to make two rather than create a bunch of methods etc...
    max_bound = 0 #this is the highest value in the list 'bound' which was passed
    for i in range(len(bound)):
        if bound[i] >= max_bound:
            max_bound = bound[i]
                
    for k in range(max_bound+1): #+1 to get around pythons counting
        ranges = [] #empty it
        for i in range(len(bound)): #increases the bound on the range on each entry, this means we dont check silly things like (0,0,...,200) 
            if k <= bound[i]:
                ranges.append(range(0,k+1))
            else:
                ranges.append(range(0,bound[i]+1)) #+1 for python counting, spent about an hour trying to debug this...
        
        for xx in itertools.product(*ranges): #test every xx in our list ranges that the "for i" loop just created for us
            successes = 0 #define a counter which resets for each xx
            for t in range(len(RHS)): #we will check each set of coeffs at 'length(RHS)' sets of evals of the schur functions. Creating a nxn linear system (in general)
                r1 = round(random.uniform(0.3,2),3)
                r2 = 0
                r3 = 0
                while r2 == 0: #if the random numbers are not distinct it all breaks down since we divide by det above in our definition, so need to do this
                    hold2 = round(random.uniform(0.3,2),3)
                    if hold2 != r1:
                        r2 = hold2
                while r3 == 0:
                    hold3 = round(random.uniform(0.3,2),3)
                    if hold3 != r1 and hold3 != r2:
                        r3 = hold3
                        
                product = 1 #this will be our varaible for the value of s(aa).s(bb) at above values. It resets at each iterate.
                valuelist = [] #this will a list of the evaluations of the schur polynomials of our possible terms evaluated at random values above. Resets each iterate.
                for i in range(len(LHS)):
                    product *= s(LHS[i][0],LHS[i][1],LHS[i][2],1,r1,r2,r3)
                for j in range(len(RHS)):
                    valuelist.append(s(RHS[j][0],RHS[j][1],RHS[j][2],1,r1,r2,r3))
                if coeff_checker(product,valuelist,xx) == True: #if it works for the random set of values, add 1 to successes and go again.
                    successes += 1
                if coeff_checker(product,valuelist,xx) == False: #if it doesnt work for a set of rand vals we can scrap the coeffs as they dont work, no need to keep checking.
                    break
            if successes == len(RHS): #if coeffs worked for every set of random numbers we tested it at, that is length(RHS) of them then return those coeffs.
                return xx
    return "DIDNT WORK, DIDNT WORK, DIDNT WORK" #very bad news...


def LittleRich(aa,bb): 
    sum_of_vals = aa[0]+aa[1]+aa[2]+bb[0]+bb[1]+bb[2]
    maxtop = aa[0] + bb[0]
    allWterms = []  #these will be all the W(a,b,c) reps with a,b,c<maxtop, a very weak restriction to start but it's something, filled with for loop right below.      
    possWterms =[]  #these will be W(a,b,c) terms from allWterms that pass some stricter conditions.
    for i in range(3): #this works for 3d weights only.               
        allWterms.append(range(0,maxtop+1)) ##+1 because range isnt an inclusive function
    
    for xs in itertools.product(*allWterms): #a variety of conditions to limit potential W terms in the linear combo.
            if xs[0]>=xs[1]>=xs[2] and xs[0]>=aa[0] and xs[0]>=bb[0] and xs[1]>=aa[1] and xs[1]>=bb[1] and xs[2] >= aa[2] and xs[2]>=bb[2] and xs[0]+xs[1]+xs[2]==sum_of_vals:
                possWterms.append(xs)
    #with any luck the list of possible W terms doesnt contain that many terms that dont appear, obviously I cant guarantee this.
    #there is a (very loose) bound on the coeffs of each one, namely dimW(aa)dimW(bb)/dim(xs[i])
    coeff_bound =[]
    for i in range(len(possWterms)):
        coeff_bound.append(int(floor(W_dim(aa)*W_dim(bb)/W_dim(possWterms[i]))))  
    
    #coeff_finder above finds the coeffs and ensures they are correct by testing potential correct ones at multiple points, creating a 'nxn linear system' (well 99.99% of the time)
    coeffs_with_0s = []
    coeffs_with_0s = LR_coeff_finder([aa,bb],possWterms,coeff_bound) #this is a list that most likely has 0s in it. We'll remove them now.
    actualWterms =[]
    actualcoeffs = []
    for i in range(len(coeffs_with_0s)): #removes terms with coeff zero, they dont appear in the expansion.
        if coeffs_with_0s[i] != 0:
            actualWterms.append(list(possWterms[i]))
            actualcoeffs.append(coeffs_with_0s[i])

    if len(actualWterms) - len(actualcoeffs) != 0: #just a small check right at the end, this should ALWAYS BE ZERO, if not we have huge problems.
        return "Error: number of coeffs doesnt equal number of terms"
    
    return [actualcoeffs,actualWterms] #returns a list of 2 lists, first list is of positive integers (multiplicity), second is of weights [a,b,c]    

def options_menu():
    print("")
    print("What would you like to do:")
    print("")
    print("1) Decompose W into a linear combination of (not necessarily irreducible) Fs. (ie perform Strong Linkage and find the coefficients)")
    print("") 
    print("2) Perform Steinberg's theorem on a F rep (ie write an F rep as a tensor product of F irreps)")
    print("")
    print("3) Perform the Littlewood-Richardson rule on two W reps")
    print("")
    print("4) Write W as a linear combination of irreducible representations (only works for W in p^2-restricted region)")
    print("")
    print("5) Write FxF as a linear combination of irreducible representations (only works for both weights in the p-restricted region)")
    print("")
    print("Please input one of the numbers, depending on what you wish to do.")
    print("(Type 'end' at any stage to, reasonably enough, end the program!)")
    print("")
    input1 = input("Input : ")
    print("")
    print("-----------------------------------------------------------------------------------------------------------------------------------------")
    return input1

def retry():
    print("")
    print("Sorry I didn't understand that, returning to options menu")
    print("")
    print("-----------------------------------------------------------------------------------------------------------------------------------------")
    main()

def thankyou():
    print("")
    print("Thank you for using the Toolkit! Goodbye :)")
    print("")
    sys.exit()
    
def finalask():
    print("")
    print("Would you like to go back to the menu or end? (type 'menu' or 'end')")    
    print("")
    input1 = input("Input : ")
    print("")
    print("-----------------------------------------------------------------------------------------------------------------------------------------")
    if [char.lower() for char in input1 if char.isalpha()] == [c.lower() for c in "menu" if c.isalpha()]: #not case sensitive and doesnt care about non-alphabet
        main()
    if [c.lower() for c in input1 if c.isalpha()] == [c.lower() for c in "end" if c.isalpha()]:
        thankyou()
    else:
        print("")
        print("Sorry I didn't understand that, please type 'menu' to go back to the menu or 'end' to end the program.")
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")    
        finalask()

def option1(weight,prime):    
    PossibleF = SLP(weight,prime) #returns list of possible weights, theres no coeffs for this yet and some could even be zero
    PostStein=[]
    for i in range(len(PossibleF)): #for each F in PossibleF we perform Steinbergs on it and add it to a new list.
        PostStein.append(Steinberg(PossibleF[i],prime)) 

    PostF_to_W=[]
    for i in range(len(PostStein)): #for every F in the PostStein nested list this changes it to combo of W terms, prepping for using schur polys to find coeffs
        hold_list =[]
        for j in range(len(PostStein[i])):
             hold_list.append(F_to_W(PostStein[i][j],prime))
        PostF_to_W.append(hold_list)

    #note, now PostF_to_W is now a triple nested list (quad technically: if you count the weights themselves as a list)
    Fcoeffs = SLP_coeff_finder(weight,PostF_to_W,prime) #returns the coeffs, could be zeros in here.

    ans_coeffs =[] #the coeffs for the answer
    ans_Fterms =[] #the terms for the answer
    for i in range(len(Fcoeffs)):
        if Fcoeffs[i] != 0:
            ans_coeffs.append(Fcoeffs[i])
            ans_Fterms.append(PossibleF[i])

    #I kept getting random "nonetype is not subscriptable" exception which would not appear upon retrying so this is the work around.
    printline = ""
    try: 
        printline += "W(" + str(weight)+") = "
        for i in range(len(ans_Fterms)-1):
            printline += str(ans_coeffs[i])+"F("+str(ans_Fterms[i])+") + " #concating coeff and term paired up as a string
        printline += str(ans_coeffs[len(ans_Fterms)-1])+"F("+str(ans_Fterms[len(ans_Fterms)-1])+")" #just concats last coeff and term pair without adding the rhs tensor sign
    except:
        printline = option1(weight,prime)
    return printline #returns string ready to print out.                   


def option2(weight,prime):        
    
    stein = Steinberg(weight,prime)
    
    #I kept getting random "nonetype is not subscriptable" eexception which would not appear upon retrying so this is the work around. It's in each option
    printline = ""
    try: 
        printline += "F(" + str(weight)+") = "
        for i in range(len(stein)-1):
            printline += "F("+str(stein[i])+")^("+str(prime**(len(stein)-1-i))+") "+u'\u2a02'+" " #concats terms and their powers.
        printline += "F("+str(stein[len(stein)-1])+")" #just concats last term
    except:
        printline = option2(weight,prime)
    return printline #returns string ready to print out.                       

def option3(weight1,weight2):
    
    LR = LittleRich(weight1,weight2)

    #I kept getting a random "nonetype is not subscriptable" exception which would not re-appear upon retrying so this is the work around. It's in each option
    printline = ""
    try: 
        printline += "W(" + str(weight1)+")" + u'\u2a02'+" W("+str(weight2)+") = "
        for i in range(len(LR[0])-1):
           printline += str(LR[0][i])+"W("+str(LR[1][i])+") + " #concats terms and their powers.
        printline += str(LR[0][len(LR[0])-1])+"W("+str(LR[1][len(LR[0])-1])+")" #just concats last term
    except:
        printline = option3(weight1,weight2)
    return printline #returns string ready to print out.

def option4(weight,prime,raw=0): ##option 4 is laid out differently to allow for recursion.          
    if weight[0] - weight[1] < 0 or weight[0]-weight[1]>= prime**2 or weight[1]-weight[2]<0 or weight[1]-weight[2]>=prime**2:            
        print("")
        print(">>>>>You have given a weight which isn't in the p^2 restricted region. Going back to menu! To try again please press 4.<<<<<")
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")            
        main()
    if 0<=weight[0]-weight[1]<prime and 0<=weight[1]-weight[2]<prime: #p-restricted, may aswell use the short code rather than run through it all if this is the case.
        printline = ""
        some_ones = []
        for i in range(len(W_to_F(weight,prime))):
            some_ones.append(1)
        try:
            printline += "W(" + str(weight)+") = "
            for i in range(len(some_ones)-1):
                printline += str(some_ones[i])+"F("+str(W_to_F(weight,prime)[i])+") + " #concating coeff and term paired up as a string
            printline += str(some_ones[len(some_ones)-1])+"F("+str(W_to_F(weight,prime)[len(some_ones)-1])+")" #just adds last term without adding the rhs + sign            
        except:
            printline = option4(weight,prime)
        return printline 

    PossibleF = SLP(weight,prime) #returns list of possible weights, theres no coeffs for this yet and some could even be zero
    PostStein=[]
    for i in range(len(PossibleF)):
        PostStein.append(Steinberg(PossibleF[i],prime)) #for each F in PossibleF we perform Steinbergs on it and add it to a new list
    PostF_to_W=[]
    for i in range(len(PostStein)): #for every F in the PostStein nested list this changes it to combo of W terms, prepping for using schur polys to find coeffs
        hold_list =[]
        for j in range(len(PostStein[i])):
             hold_list.append(F_to_W(PostStein[i][j],prime))
        PostF_to_W.append(hold_list)
    #note, now PostF_to_W is now a triple nested list (quad technically: if you count the weights themselves as a list)
    Fcoeffs = SLP_coeff_finder(weight,PostF_to_W,prime) #returns the coeffs, could be zeros in here.
    good_coeffs =[]
    good_Fterms = []
    for i in range(len(Fcoeffs)):
        if Fcoeffs[i] != 0:
            good_coeffs.append(Fcoeffs[i])
            good_Fterms.append(PostF_to_W[i]) #return list, first entry are the inputted weights, second entry is the answer.
    #we want to expand/rearrange this for loop a bit before doing LR, so let's do that.
    #this is why it only works for p^2, it got too confusing to work with >2 terms in each good_Fterms[i] bc i couldnt just use a 'j' and 'k' for loop.
    
    expansion = []
    expansion_coeffs = []
    for i in range(len(good_Fterms)): #we want to run through every "good term"
        if len(good_Fterms[i]) ==1:#this implies theres only the term in the Steinberg decomp, so no LR needed.
            for j in range(len(good_Fterms[i][0])): #we want to go into this term (it might be in the form W(wgt) - W(wgt)
                expansion.append([good_Fterms[i][0][j]]) #add it to our new list, this new list will be double nested hence the extra [ ](triple if you count weights as list).
                expansion_coeffs.append(good_coeffs[i]*(-1)**j) #add coeff to list, the (-1)**j now accounts for whether the W term had a minus at the front.
        if len(good_Fterms[i]) == 2:
            for j in range(len(good_Fterms[i][0])): #length of first term in the list, this term may be of form 'W-W',same with second term in list just below
                for k in range(len(good_Fterms[i][1])): #length of second term in the list, this first and second term get tensored together.
                    expansion.append([good_Fterms[i][0][j],good_Fterms[i][1][k]]) #same as above, i need to manually add an extra [ ] 
                    expansion_coeffs.append(good_coeffs[i]*(-1)**(j+k)) #same idea as above, this time the j+k brings into play the minus signs.
        if len(good_Fterms[i]) > 2:
            print("shouldnt be in here as this would imply 3 tensor decomp which I havent coded for. Running will cease. (im bad!)")
            sys.exit()
  
    #now we have two new lists, expansion which looks like [ [weight,weight],[weight,weight],[weight],...] and expansion_coeffs which is just a list of coeffs.
    #if length of any term in expansion is 2 we need to perform LR, this also involves adding more coeffs (the same one multi times) so let's create a further two new lists.
    almost_terms =[] #called so because we can see the light at the end of the tunnel. almost there!
    almost_coeffs =[]
    for i in range(len(expansion)):
        if len(expansion[i]) == 1:
            almost_terms.append(expansion[i][0]) #add the only weight in expansion[i1 to new list, new list is now a list of weights. no more nesting :D
            almost_coeffs.append(expansion_coeffs[i]) #add the coeff too one time, nothing crazy happening here
        if len(expansion[i]) == 2: #we need to use LR at last
            result = LittleRich(expansion[i][0],expansion[i][1]) #remember this returns a list of lists, 1st list the coeffs 2nd the terms.
            for j in range(len(result[1])): #note, could write result[0] here also, they have the same length.
                almost_terms.append(result[1][j]) #run through appending all terms of expansion to the big list.
                almost_coeffs.append(result[0][j]*expansion_coeffs[i]) #here append the product LR coeff and the coeff from the SLP coeff that carries through.
    #for W terms outside p-restricted region is use recursion.
    
    #in this block we wish to now convert this list of Ws into irrep Fs, this will use alcove position and maybe even recursion
    master_terms = []
    master_coeffs =[]
    for i in range(len(almost_terms)):
        a,b,c = almost_terms[i][0],almost_terms[i][1], almost_terms[i][2] #weights in shorthand
        if 0 <= a-b and 0 <= b-c and a-c <=prime-2 or a-b == prime - 1 and b-c <= prime - 1 or b-c == prime -1 and a-b <= prime -1: #lower alcove and boundary
            master_terms.append(almost_terms[i])
            master_coeffs.append(almost_coeffs[i])
        elif a-b < prime - 1 and b - c < prime - 1 and a-c > prime - 2:
            master_terms.append(W_to_F(almost_terms[i],prime)[0]) #append first term of the decomp of upper alcove W to master list
            master_coeffs.append(almost_coeffs[i]) #add coeff associated to W term, the coeff of both terms in decomp is 1 so no scaling/minus to worry about.
            master_terms.append(W_to_F(almost_terms[i],prime)[1]) #do the same for the second term in expansion
            master_coeffs.append(almost_coeffs[i]) #add coeff again, it's the same one as two lines above.
        else: # W is outside p-restricted, use recursion here.
            mess = option4(almost_terms[i],prime,1) #this returns decomp of the W term not in p-restricted.
            for j in range(len(mess[1])): #mess[2] has same value,could use it instead
                master_terms.append(mess[1][j]) #add each term in decomp of 'almost_terms[i]' to master list
                master_coeffs.append(mess[0][j]*almost_coeffs[i]) #all coeffs in this decomp are scaled by almost_coeff[i]. (think substutition, we replaced the bad W term with mess, a list of good ones,scaler still acts)
            
    #master terms and coeffs is technically correct, altho it contains minus terms (of course, there is a corres. +ve to cancel it) and repeated terms, so lets tidy the lists up
    tidy_terms =[]
    tidy_coeffs=[]
    for i in range(len(master_terms)):
        multiplicity = 0 #to be reset every loop
        if len(tidy_terms) == 0: #empty
            for j in range(len(master_terms)):
                if (master_terms[i][0] - master_terms[j][0] == master_terms[i][1] - master_terms[j][1] #first says equal up to spin. fourth => iso.
                      == master_terms[i][2] - master_terms[j][2] and (master_terms[i][0] - master_terms[j][0]) % (prime-1) == 0):
                    #you get into this if statement iff master_terms[i] is isomorphic to master_terms[j]
                    multiplicity += master_coeffs[j] #if master_term[i] = master_term[j] then add term[j]'s multiplicity to the overall value
            if(multiplicity != 0): #if the multiplicity is zero (ie some terms cancelled) then its silly to add it (altho not incorrect)
                    reduced_term = master_terms[i].copy() #set term to be added equal to master_terms[i]
                    while reduced_term[0] >= prime -1 and reduced_term[1] >= prime -1 and reduced_term[2] >= prime -1: #if true then reduce term via isomorphism
                        reduced_term[0] -= prime -1
                        reduced_term[1] -= prime -1
                        reduced_term[2] -= prime -1
                    tidy_terms.append(reduced_term) #now add term to tidy_terms along with it's non zero multiplicity 
                    tidy_coeffs.append(multiplicity)                    
        else: #tidy_terms is not empty, we must be careful not to repeat terms.
            notiso = 0
            for u in range(len(tidy_terms)):                    
                if ( master_terms[i][0] - tidy_terms[u][0] != master_terms[i][1] - tidy_terms[u][1] or #conditions for not iso
                        master_terms[i][1] - tidy_terms[u][1] != master_terms[i][2] - tidy_terms[u][2] or
                        (master_terms[i][0] - tidy_terms[u][0]) % (prime-1) !=0):
                    notiso += 1
            if notiso == len(tidy_terms): #happens if master_terms[i] is not iso to any of the terms in tidy_terms, so let's add it in alond with multiplicty.
                for j in range(len(master_terms)):
                    if (master_terms[i][0] - master_terms[j][0] == master_terms[i][1] - master_terms[j][1] #first says equal up to spin. fourth => iso.
                          == master_terms[i][2] - master_terms[j][2] and (master_terms[i][0] - master_terms[j][0]) % (prime-1) == 0):
                        #you get into this if statement iff master_terms[i] is isomorphic to master_terms[j]
                        multiplicity += master_coeffs[j]
                if(multiplicity != 0): #if the multiplicity is zero (ie some terms cancelled) then its silly to add it (altho not incorrect)
                    reduced_term = master_terms[i].copy() #set term to be added equal to master_terms[i]
                    while reduced_term[0] >= prime -1 and reduced_term[1] >= prime -1 and reduced_term[2] >= prime -1: #if true then reduce term via isomorphism
                        reduced_term[0] -= prime -1
                        reduced_term[1] -= prime -1
                        reduced_term[2] -= prime -1
                    tidy_terms.append(reduced_term) #now add term to tidy_terms along with it's non zero multiplicity 
                    tidy_coeffs.append(multiplicity)
        
    #I kept getting random "nonetype is not subscriptable" exception which would not appear upon retrying so this is the work around.
    if raw == 0: #for recursion I want to be able to call option(4) again and get the 'raw' output ie the lists, not a string.
        printline = ""
        try:
            printline += "W(" + str(weight)+") = "
            for i in range(len(tidy_terms)-1):
                printline += str(tidy_coeffs[i])+"F("+str(tidy_terms[i])+") + " #concating coeff and term paired up as a string
            printline += str(tidy_coeffs[len(tidy_coeffs)-1])+"F("+str(tidy_terms[len(tidy_terms)-1])+")" #just adds last term without adding the rhs + sign            
        except:
            printline = option4(weight,prime)
        return printline
    if raw == 1:
        return [tidy_coeffs,tidy_terms]

def option5(weight1,weight2,prime):
        PostF_to_W = [] #this will be a list of W weights depending on whether the Fs lie in lower or upper alcove
        PostF_to_W.append(F_to_W(weight1,prime))
        PostF_to_W.append(F_to_W(weight2,prime))
        
        expansion=[] #a list of terms after using districution/linearity properties of tensor on (W-w)x(W-W)
        expansion_sign = [] #to keep track of the +/- signs.
        if len(PostF_to_W[0]) == 1:
            for i in range(len(PostF_to_W[1])):
                expansion.append([PostF_to_W[0][0],PostF_to_W[1][i]])
                expansion_sign.append((-1)**i) #add a +/- depending on whether the W term in the second tensor has a plus or minus infront of it.
        if len(PostF_to_W[0]) == 2:
            for i in range(len(PostF_to_W[0])):                
                for j in range(len(PostF_to_W[1])):
                    expansion.append([PostF_to_W[0][i],PostF_to_W[1][j]])
                    expansion_sign.append((-1)**(i+j)) #add a +/- depending on sign of WxW after expanding (W-W)x(W-W) with linearity.
        final_terms = []
        final_coeffs =[]
        for i in range(len(expansion)):            
            coeff_term_lists = LittleRich(expansion[i][0],expansion[i][1])
            for j in range(len(coeff_term_lists[1])): #runs through every term in the result of LR to see what to do with it.
                if coeff_term_lists[1][j][0] - coeff_term_lists[1][j][1] >= prime or coeff_term_lists[1][j][1] - coeff_term_lists[1][j][2] >= prime: #if outside p-restricted call option4
                    term_decomp = option4(coeff_term_lists[1][j],prime,1) #decomp W term into irrep F terms by calling option4().
                    for k in range(len(term_decomp[1])):
                        final_terms.append(term_decomp[1][k]) #add each irrep F coming from option4() being called to the final list.
                        final_coeffs.append(expansion_sign[i]*coeff_term_lists[0][j]*term_decomp[0][k]) #add coeff also, remembering the =/- and coeff from LR
                #if the above isnt true were in the p-restricted region, now to decided whether were in upper or low alcove
                elif (coeff_term_lists[1][j][0] - coeff_term_lists[1][j][2] > prime - 2 and
                          coeff_term_lists[1][j][0] - coeff_term_lists[1][j][1] < prime - 1 and
                          coeff_term_lists[1][j][1] - coeff_term_lists[1][j][2] < prime - 1): #upper alcove conds
                    term_decomp = W_to_F(coeff_term_lists[1][j],prime)
                    for k in range(2):
                        final_terms.append(term_decomp[k]) #add each irrep F coming from option4() being called to the final list.
                        final_coeffs.append(expansion_sign[i]*coeff_term_lists[0][j])
                else: #we must be in the lower alcove or boundary.
                    final_terms.append(coeff_term_lists[1][j]) #add the single term since W = F in lower alcove
                    final_coeffs.append(expansion_sign[i]*coeff_term_lists[0][j])

        #here we want to group terms that are isomorphic and have them all together with a combined coeff.            
        tidy_terms =[]
        tidy_coeffs =[]
        for i in range(len(final_terms)):
            multiplicity = 0
            if len(tidy_terms) == 0: #empty
                for j in range(len(final_terms)):
                    if (final_terms[i][0] - final_terms[j][0] == final_terms[i][1] - final_terms[j][1] #first says equal up to spin. second cond => iso.
                          == final_terms[i][2] - final_terms[j][2] and (final_terms[i][0] - final_terms[j][0]) % (prime-1) == 0):
                        #you get into this if statement iff final_terms[i] is isomorphic to final_terms[j]
                        multiplicity += final_coeffs[j] #if final_terms[i] = final_terms[j] then add term[j]'s multiplicity to the overall value
                if(multiplicity != 0): #if the multiplicity is zero (ie some terms cancelled) then its silly to add it (altho not incorrect)
                        reduced_term = final_terms[i].copy() #set term to be added equal to master_terms[i]
                        while reduced_term[0] >= prime -1 and reduced_term[1] >= prime -1 and reduced_term[2] >= prime -1: #if true then reduce term via isomorphism
                            reduced_term[0] -= prime -1
                            reduced_term[1] -= prime -1
                            reduced_term[2] -= prime -1
                        tidy_terms.append(reduced_term) #now add term to tidy_terms along with it's non zero multiplicity 
                        tidy_coeffs.append(multiplicity)
            else: #tidy_terms is not empty, we must be careful not to repeat terms.
                notiso = 0
                for u in range(len(tidy_terms)):                    
                    if ( final_terms[i][0] - tidy_terms[u][0] != final_terms[i][1] - tidy_terms[u][1] or #conditions for not iso
                            final_terms[i][1] - tidy_terms[u][1] != final_terms[i][2] - tidy_terms[u][2] or
                            (final_terms[i][0] - tidy_terms[u][0]) % (prime-1) !=0):
                        notiso += 1
                if notiso == len(tidy_terms): #happens if final_terms[i] is not iso to any of the terms in tidy_terms, so let's add it in along with multiplicty.
                    for j in range(len(final_terms)):
                        if (final_terms[i][0] - final_terms[j][0] == final_terms[i][1] - final_terms[j][1] #first says equal up to spin. fourth => iso.
                              == final_terms[i][2] - final_terms[j][2] and (final_terms[i][0] - final_terms[j][0]) % (prime-1) == 0):
                            #you get into this if statement iff master_terms[i] is isomorphic to master_terms[j]
                            multiplicity += final_coeffs[j]
                    if(multiplicity != 0): #if the multiplicity is zero (ie some terms cancelled) then its silly to add it (altho not incorrect)
                            reduced_term = final_terms[i].copy() #set term to be added equal to master_terms[i]
                            while reduced_term[0] >= prime -1 and reduced_term[1] >= prime -1 and reduced_term[2] >= prime -1: #if true then reduce term via isomorphism
                                reduced_term[0] -= prime -1
                                reduced_term[1] -= prime -1
                                reduced_term[2] -= prime -1
                            tidy_terms.append(reduced_term) #now add term to tidy_terms along with it's non zero multiplicity 
                            tidy_coeffs.append(multiplicity)           

        printline =""
        try:
            printline += "F(" + str(weight1)+")"+u'\u2a02'+"F("+str(weight2)+" = "
            for i in range(len(tidy_terms)-1):
                printline += str(tidy_coeffs[i])+"F("+str(tidy_terms[i])+") + " #concating coeff and term paired up as a string
            printline += str(tidy_coeffs[len(tidy_coeffs)-1])+"F("+str(tidy_terms[len(tidy_terms)-1])+")" #just adds last term without adding the rhs + sign            
        except:
            printline = option5(weight1,weight2,prime)
        return printline
                                     
            
def main():
    input_options = options_menu()
    
    if input_options == "1":
        print("")    
        print("**********************************************************************************************************")
        print("O P T I O N  1: W(a,b,c) D E C O M P O S I T I O N  I N T O  (N O T  N E C E S S A R I L Y  I R R E P S) Fs ")
        print("**********************************************************************************************************")
        print("")
        print("Please input the rep W(a,b,c) you wish to decompose and the prime in the form 'a,b,c,p' or type 'back' to go back to the options")
        print("")
        print("Input : ",end = ' ')
        input_list = list(tuple(x.strip() for x in input().split(','))) #takes a list of inputs where we seperate elemnts via ',' in the input
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")
        if input_list[0] == 'back': #if user typed back in would be the first entry of the list, send them back 
            main()
        elif input_list[0] == 'end': #ditto above but for end
            thankyou()
        elif len(input_list) == 4: 
            try:
                weight = [int(input_list[0]),int(input_list[1]),int(input_list[2])] #are the inputs ints?? this try tells us the answer.
                prime = int(input_list[3])
            except:#not all ints, make user try again
                retry()
        else:#this covers inputs which dont have 4 entries, this could be a variety of things so best to just make them try again rather than try guess what they meant.
            retry()
        if weight[0] < weight[1] or weight[1]<weight[2]:
            print("")            
            print("You have given a weight which is not in the dominant region, returning to menu, please try again.")
            print("")
            print("-----------------------------------------------------------------------------------------------------------------------------------------")
            main()            
        #if all the above are ok we'll go for an answer
        ans = option1(weight,prime) 
        print("")    
        print("**************")
        print(" A N S W E R ")
        print("**************")
        print("")
        print(ans)
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")
                   
    elif input_options == "2":
        print("")    
        print("*********************************************************************")
        print(" O P T I O N  2 :  S T E I N B E R G ' S  T E N S O R  T H E O R E M ")
        print("*********************************************************************")
        print("")
        print("Please input the weight of the F(a,b,c) rep you wish to decompose and the prime in the form 'a,b,c,p' or type 'back' to go back to the options")
        print("")
        print("Input : ",end = ' ')
        input_list = list(tuple(x.strip() for x in input().split(','))) #takes a list of inputs where we seperate elemnts via ',' in the input
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")
        if input_list[0] == 'back': #if user typed back in would be the first entry of the list, send them back 
            main()
        elif input_list[0] == 'end': #ditto above but for end
           thankyou()
        elif len(input_list) == 4: 
            try:
                weight = [int(input_list[0]),int(input_list[1]),int(input_list[2])] #are the inputs ints?? this try tells us the answer.
                prime = int(input_list[3])
            except:#not all ints, make user try again
                retry()
        else:#this covers inputs which dont have 4 entries, this could be a variety of things so best to just make them try again rather than try guess what they meant.
            retry()
        if weight[0] < weight[1] or weight[1]<weight[2]:
            print("")            
            print("You have given a weight which is not in the dominant region, returning to menu, please try again.")
            print("")
            print("-----------------------------------------------------------------------------------------------------------------------------------------")
            main()

        #if we get pass all of the above its time to find an answer and print it
        ans = option2(weight,prime)
        print("")    
        print("**************")
        print(" A N S W E R ")
        print("**************")
        print("")
        print(ans)
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")
                
    elif input_options == "3": 
        print("")    
        print("***************************************************************************************")
        print(" O P T I O N  3 :  L I T T L E W O O D - R I C H A R D S O N  O N  W(a,b,c) "+u'\u2a02'+" W(d,e,f)")
        print("***************************************************************************************")
        print("")
        print("Please input the weights of two reps W(a,b,c) and W(d,e,f) in the form 'a,b,c,d,e,f' or type 'back' to go back to the options")
        print("")
        print("Input : ",end = ' ')
        input_list = list(tuple(x.strip() for x in input().split(','))) #takes a list of inputs where we seperate elemnts via ',' in the input
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")
        if input_list[0] == 'back': #if user typed back in would be the first entry of the list, send them back 
            main()
        elif input_list[0] == 'end': #ditto above but for end
           thankyou()
        elif len(input_list) == 6: 
            try:
                weight1 = [int(input_list[0]),int(input_list[1]),int(input_list[2])]
                weight2 = [int(input_list[3]),int(input_list[4]),int(input_list[5])] #are the inputs ints?? this try tells us the answer.            
            except:#not all ints, make user try again
                retry() 
        else:
            retry()            
        if weight1[0] < weight1[1] or weight1[1]<weight1[2] or weight2[0] < weight2[1] or weight2[1]<weight2[2]:
            print("")
            print("You have entered a weight which isn't in the dominant region. Going back to menu.")
            print("")
            print("-----------------------------------------------------------------------------------------------------------------------------------------")
            main()
        if weight1[0] < 0 or weight1[1] < 0 or weight1[2]< 0 or weight2[0] < 0 or weight2[1] < 0 or weight2[2]< 0:
            print("")
            print("Invalid input, please ensure all entries are >= 0. Going back to menu.")
            print("")
            print("-----------------------------------------------------------------------------------------------------------------------------------------")
            main()

        ans = option3(weight1,weight2)
        print("")    
        print("**************")
        print(" A N S W E R ")
        print("**************")
        print("")
        print(ans)
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")
                
    elif input_options == "4": #done differently because i wanted option4(paras) to be available for recursion
        print("")
        print("****************************************************************************************************")
        print("O P T I O N  4 : F U L L  D E C O M P O S I T I O N  O F  W(a,b,c) I N T O  S U M  O F  I R R E P S")
        print("****************************************************************************************************")
        print("")
        print("Input the weight of a rep W(a,b,c) in the p^2 restricted region you wish to decompose and the prime in the form 'a,b,c,p' or type 'back' to go back to the options")
        print("")
        print("Input : ",end = ' ')
        input_list = list(tuple(x.strip() for x in input().split(','))) #takes a list of inputs where we seperate elemnts via ',' in the input
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")            
        if input_list[0] == 'back': #if user typed back in would be the first entry of the list, send them back 
            main()
        elif input_list[0] == 'end': #ditto above but for end
           thankyou()
        elif len(input_list) == 4: 
            try:
                weight = [int(input_list[0]),int(input_list[1]),int(input_list[2])]
                prime = int(input_list[3]) #are the inputs ints?? this try tells us the answer.            
            except:#not all ints, make user try againprint("")
                retry()
        else:#this covers inputs which dont have 6 entries, this could be a variety of things so best to just make them try again rather than try guess what they meant.
            retry()
        if weight[0] < weight[1] or weight[1]<weight[2] or weight[0]-weight[1] >= prime**2 or weight[1]-weight[2] >= prime**2:
            print("")            
            print("You have given a weight which is not in the p^2-restricted region, returning to menu, please try again.")
            print("")
            print("-----------------------------------------------------------------------------------------------------------------------------------------")
            main()
        ans = option4(weight,prime,0) #the 0 means we'll get a nicely laid out string and not a list of lists.
        print("")    
        print("**************")
        print(" A N S W E R   ")
        print("**************")        
        print("")        
        print(ans)
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")

    elif input_options == "5": 
        print("")
        print("*****************************************************************************************")
        print("O P T I O N  5 : D E C O M P O S I T I O N  O F  F(a,b,c)"+u'\u2a02'+" F(d,e,f) I N T O  I R R E P S")
        print("*****************************************************************************************")
        print("")
        print("Input the weights of two reps F(a,b,c), F(d,e,f) in the p - restricted region you wish to decompose and the prime in the form 'a,b,c,d,e,f,p' or type 'back' to go back to the options")
        print("")
        print("Input : ",end = ' ')
        input_list = list(tuple(x.strip() for x in input().split(','))) #takes a list of inputs where we seperate elemnts via ',' in the input
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")            
        if input_list[0] == 'back': #if user typed back in would be the first entry of the list, send them back 
            main()
        elif input_list[0] == 'end': #ditto above but for end
           thankyou()
        elif len(input_list) == 7: 
            try: #are the inputs ints?? this try catch ensures they are before going ahead. 
                weight1 = [int(input_list[0]),int(input_list[1]),int(input_list[2])]
                weight2 = [int(input_list[3]),int(input_list[4]),int(input_list[5])]
                prime = int(input_list[6])            
            except:#not all ints, make user try againprint("")
                retry()
        else:#this covers inputs which dont have 7 entries, this could be a variety of things so best to just make them try again rather than try guess what they meant.
            retry()
        if (weight1[0] < weight1[1] or weight1[1] < weight1[2] or weight1[0]-weight1[1] >= prime or weight1[1]-weight1[2] >= prime or
            weight2[0] < weight2[1] or weight2[1] < weight2[2] or weight2[0]-weight2[1] >= prime or weight2[1]-weight2[2] >= prime):
            print("")            
            print("You have given a weight which is not in the p-restricted region, returning to menu, please try again.")
            print("")
            print("-----------------------------------------------------------------------------------------------------------------------------------------")
            main()
        ans = option5(weight1,weight2,prime) #the 0 means we'll get a nicely laid out string and not a list of lists.
        print("")    
        print("**************")
        print(" A N S W E R   ")
        print("**************")        
        print("")        
        print(ans)
        print("")
        print("-----------------------------------------------------------------------------------------------------------------------------------------")
        
    elif [c.lower() for c in input_options if c.isalpha()] == [c.lower() for c in "end" if c.isalpha()]:
        thankyou()
    else: #didnt recognise input
        retry()
        main()      
    #once you finish a cal this is asked
    finalask()
    sys.exit() 

print("-----------------------------------------------------------------------------------------------------------------------------------------")
print("")
print("*****************************************************************************************")
print("*  W E L C O M E  T O  T H E  G L 3(F p)  R E P R E S E N T A T I O N S  T O O L K I T  *")
print("*****************************************************************************************")
main()
