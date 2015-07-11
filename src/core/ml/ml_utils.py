from sklearn.externals import joblib

def SaveModel(model, str_of_md):
    joblib.dump(model, str_of_md)
