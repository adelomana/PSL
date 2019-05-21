import numpy

def randomCohortMaker():

    cohort=[]
    for i in range(trialSize):
        p=numpy.random.random()
        if p <= responderLikelihood:
            cohort.append(True)
        else:
            cohort.append(False)

    return cohort


# MAIN
responderLikelihood=32/870
trialSize=25

rol=1e9
requirement=10

trial=0
successes=0
while trial < rol:
    trial=trial+1
    cohort=randomCohortMaker()

    if trial % (rol/100) == 0:
        print('\t working on trial {:.2E}'.format(trial))

    if sum(cohort) >= requirement:
        successes=successes+1
        frequency=successes/trial
        print('FOUND: trial {}; successes {}; frequency {:.2E}'.format(trial,successes,frequency))
print(successes)
