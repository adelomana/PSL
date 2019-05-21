import numpy

# MAIN
likelihoodA=32/870
likelihoodB=1-likelihoodA
cohortSize=25
requestedSize=10

iterations=1
draws=['A','B']
while iterations < cohortSize:
    iterations=iterations+1

    new=[]
    for draw in draws:
        new.append(draw+'A')
        new.append(draw+'B')
    draws=new
    #print(iterations,draws,len(draws))
    print('Iteration: {}; combinations: {}'.format(iterations,len(draws)))
    
    # compute the number likelihood of each draw
    total=0; requestedLikelihood=0
    for draw in draws:
        likelihood=1
        for element in draw:
            if element == 'A':
                likelihood=likelihood*likelihoodA
            if element == 'B':
                likelihood=likelihood*likelihoodB

        # taking good draws
        if draw.count('A') >= requestedSize:
            requestedLikelihood=requestedLikelihood+likelihood
            accepting=True
        else:
            accepting=False
            
        # checking full likelihood and prints    
        total=total+likelihood
        #print(draw,likelihood,accepting)
    print('total: {:.5E}'.format(total))
    print('requestedLikelihood: {:.5E}'.format(requestedLikelihood))
    print()
    

