from sklearn.externals import joblib

def SaveModel(mdl, str_of_md):
    joblib.dump(mdl, str_of_md, compress=3)

def LoadModel(str_of_md):
    return joblib.load(str_of_md)


def LassoEval(model, x, columns):
    return model.predict(x[:,columns]).ravel()

def LogitEval(model, x, columns):
    return model.predict_proba(x[:,columns])[:,1].ravel()
