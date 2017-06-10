import sys

from sklearn.externals import joblib
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier

scaler = joblib.load('scalermlp.sav')
clf = joblib.load('mlp.sav')

'''
while True:
    aa = []
    f = []
    for line in sys.stdin:
        split = line.split()
        aa.append([split[0], split[1]])
        f.append([float(x) for x in split[2:]])
    if (f == []): continue
    X = scaler.transform(f)
    pred = clf.predict_proba(X)
    for i in range(len(f)):
        print("%s %s %f" % (aa[i][0], aa[i][1], pred[i][1]))
    sys.stdout.flush()
'''

while True:
    for line in sys.stdin:
        print(clf.predict_proba(
                scaler.transform(
                    [[float(x) for x in line.split()]])
            )
              [0][1])
        sys.stdout.flush()
