import pandas as pd
import xlsxwriter
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import AdaBoostRegressor,GradientBoostingRegressor,BaggingRegressor,ExtraTreesRegressor,RandomForestRegressor 
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, Matern, WhiteKernel, Exponentiation, Product, Sum
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error, mean_squared_error
import numpy as np
import pickle as pkl

def Optimized_payan_name_Tg_nano_composite():
    RF_GB_models = [
              RandomForestRegressor(n_estimators = 217, min_samples_split = 0.093840625),
              RandomForestRegressor(n_estimators = 272, min_samples_split = 0.031227635),
              RandomForestRegressor(n_estimators = 350, min_samples_split = 0.039243403),
              RandomForestRegressor(n_estimators = 107, min_samples_split = 0.126663933),
              RandomForestRegressor(n_estimators = 140, min_samples_split = 0.024252354),
              GradientBoostingRegressor(n_estimators = 178, min_samples_split = 0.18758125 ),
              GradientBoostingRegressor(n_estimators = 286, min_samples_split = 0.049154266),
              GradientBoostingRegressor(n_estimators = 167, min_samples_split = 0.193593156),
              GradientBoostingRegressor(n_estimators = 255, min_samples_split = 0.138509697),
              GradientBoostingRegressor(n_estimators = 286, min_samples_split = 0.049154266)
              ]
    ET_models = [
              # extra trees optimized results
              ExtraTreesRegressor(n_estimators=250, min_samples_split= 0.163734361),
              ExtraTreesRegressor(n_estimators=224, min_samples_split= 0.143697764),
              ExtraTreesRegressor(n_estimators=224, min_samples_split= 0.143697764),
              ExtraTreesRegressor(n_estimators=241, min_samples_split= 0.10346333 ),
              ExtraTreesRegressor(n_estimators=224, min_samples_split= 0.143697764),
              ]
    GPR_models = [
              # gaussian process rbf optimized results
              GaussianProcessRegressor(kernel = RBF(length_scale = 2.535937332)),
              GaussianProcessRegressor(kernel = RBF(length_scale = 2.536145069)),
              GaussianProcessRegressor(kernel = RBF(length_scale = 0.88173322 )),
              GaussianProcessRegressor(kernel = RBF(length_scale = 1.39160148 )),
              GaussianProcessRegressor(kernel = RBF(length_scale = 2.536145069)),
              # Gaussian Process Exponent RBF optimized results
              GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 2.709138934),exponent = 0.245295861)),
              GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 1.510927643),exponent = 0.262880944)),
              GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 1.011911938),exponent = 0.117456255)),
              GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 2.752461192),exponent = 0.823807167)),
              GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 1.510927643),exponent = 0.262880944)),
              # Gaussian Process Matern optimized results
              GaussianProcessRegressor(kernel=Matern(length_scale = 2.548434499, nu = 0.5)),
              GaussianProcessRegressor(kernel=Matern(length_scale = 0.480590624, nu = 0.5)),
              GaussianProcessRegressor(kernel=Matern(length_scale = 0.480590624, nu = 0.5)),
              GaussianProcessRegressor(kernel=Matern(length_scale = 0.218565733, nu = 0.5)),
              GaussianProcessRegressor(kernel=Matern(length_scale = 0.480590624, nu = 0.5)),
              # Gaussian Process Exponent Matern optimized results
              GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 2.775323549, nu=  1.5), exponent = 0.246048782)),
              GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 2.512309459, nu=  2.5), exponent = 0.246429492)),
              GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 2.512309459, nu=  2.5), exponent = 0.246429492)),
              GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 1.738734569, nu=  2.5), exponent = 0.010661179)),
              GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 2.512309459, nu=  2.5), exponent = 0.246429492))
              ]
    
    filename = "dataset"
    dataset = pkl.load(open("data/"+filename, "rb"))

    X_fields = ["graphene_w","CH","C=OH","OH","NH","B","BCH","hydrogen_bond","polarity"]
    X = normalize(np.array(dataset[["P_T","graphene_w","CH","C=OH","OH","NH","B","BCH","hydrogen_bond","polarity"]]),'l2',axis=0)
    Y = np.array(dataset[["nano_composite_Tg_(K)"]])
    TEST_SIZE = 0.15
    polymer_names_with_Tg = dataset.iloc[:,[1,6]]
    
    # training rf_gb models
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=TEST_SIZE, random_state = 42)
    Xtrain = np.array(Xtrain)
    Xtest  = np.array(Xtest)
    Ytrain = np.array(Ytrain).reshape(len(Ytrain),1)
    Ytest  = np.array(Ytest).reshape(len(Ytest),1)
    model_names = [str(name).split("()")[0] for name in RF_GB_models]
    number_of_features = X.shape[1]
    number_of_models = len(model_names)
    polymer_names = ["" for _ in range(dataset.shape[0])]
    all_Y_real = np.concatenate((Ytrain,Ytest),axis=0) # for calculating MAE, MAPE and MSE
    index = 0
    for Tg in all_Y_real.flatten():
        # find the index of the polymer with the matching Tg value
        Tg_index = polymer_names_with_Tg[polymer_names_with_Tg["nano_composite_Tg_(K)"] == Tg].index
        if not Tg_index.empty:  # check if any index was found
            polymer_names[index] = polymer_names_with_Tg.iloc[Tg_index[0]]["polymer_name"]
            index += 1
        else:
            print("could not find this Tg",Tg)

    len_Xtest  = Xtest.shape[0]
    len_Xtrain = Xtrain.shape[0]
    results_train = np.zeros((len_Xtrain,number_of_models))
    results_test = np.zeros((len_Xtest,number_of_models))
    results_R2_train = np.zeros((1,number_of_models))
    results_R2_test = np.zeros((1,number_of_models))
    importances_results = np.zeros((number_of_features,number_of_models))
    MAE_results = np.zeros((1,number_of_models))
    MAPE_results = np.zeros((1,number_of_models))
    MSE_results = np.zeros((1,number_of_models))
    number_of_train_data = Ytrain.shape[0]
    
    for index, regressor in enumerate(RF_GB_models):    
        regressor.fit(Xtrain,Ytrain.ravel())
        predict_train = regressor.predict(Xtrain)
        predict_test = regressor.predict(Xtest)
        all_Y_pred = np.concatenate((predict_train,predict_test))
        MAE_results[0,index] = mean_absolute_error(all_Y_real,all_Y_pred)
        MAPE_results[0,index] = mean_absolute_percentage_error(all_Y_real,all_Y_pred)
        MSE_results[0,index] = mean_squared_error(all_Y_real,all_Y_pred)
        results_train[:,index] = regressor.predict(Xtrain) 
        results_test[:,index] = regressor.predict(Xtest) 
        r2_train = regressor.score(Xtrain,Ytrain)
        r2_test  = regressor.score(Xtest,Ytest)
        if (r2_train < 0):
            r2_train = np.abs(r2_train)
        if (r2_test < 0):
            r2_test = np.abs(r2_test)
        if (r2_train > 1):
            r2_train = 1/r2_train
        if (r2_test > 1):
            r2_test = 1/r2_test

        results_R2_train[0,index] = r2_train
        results_R2_test[0,index]  = r2_test 
        importances = permutation_importance(regressor,Xtrain,Ytrain)
        importances_results[:,index] = importances["importances_mean"]
        print(str(regressor).split("()")[0]+" complete.")
    print("saving results")
    Save_results(Ytrain,results_train,Ytest,results_test,results_R2_train,results_R2_test,MAE_results,MAPE_results,MSE_results,model_names,(polymer_names,number_of_train_data),importances_results,X_fields,"RF_GB_Optimized_predictions")
    
    # training et models
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=TEST_SIZE, random_state = 42)
    Xtrain = np.array(Xtrain)
    Xtest  = np.array(Xtest)
    Ytrain = np.array(Ytrain).reshape(len(Ytrain),1)
    Ytest  = np.array(Ytest).reshape(len(Ytest),1)
    model_names = [str(name).split("()")[0] for name in ET_models]
    number_of_features = X.shape[1]
    number_of_models = len(model_names)
    all_Y_real = np.concatenate((Ytrain,Ytest),axis=0) # for calculating MAE, MAPE and MSE
    index = 0
    for Tg in all_Y_real.flatten():
        # find the index of the polymer with the matching Tg value
        Tg_index = polymer_names_with_Tg[polymer_names_with_Tg["nano_composite_Tg_(K)"] == Tg].index
        if not Tg_index.empty:  # check if any index was found
            polymer_names[index] = polymer_names_with_Tg.iloc[Tg_index[0]]["polymer_name"]
            index += 1
        else:
            print("could not find this Tg",Tg)

    len_Xtest  = Xtest.shape[0]
    len_Xtrain = Xtrain.shape[0]
    results_train = np.zeros((len_Xtrain,number_of_models))
    results_test = np.zeros((len_Xtest,number_of_models))
    results_R2_train = np.zeros((1,number_of_models))
    results_R2_test = np.zeros((1,number_of_models))
    importances_results = np.zeros((number_of_features,number_of_models))
    MAE_results = np.zeros((1,number_of_models))
    MAPE_results = np.zeros((1,number_of_models))
    MSE_results = np.zeros((1,number_of_models))
    number_of_train_data = Ytrain.shape[0]
    
    for index, regressor in enumerate(ET_models):    
        regressor.fit(Xtrain,Ytrain.ravel())
        predict_train = regressor.predict(Xtrain)
        predict_test = regressor.predict(Xtest)
        all_Y_pred = np.concatenate((predict_train,predict_test))
        MAE_results[0,index] = mean_absolute_error(all_Y_real,all_Y_pred)
        MAPE_results[0,index] = mean_absolute_percentage_error(all_Y_real,all_Y_pred)
        MSE_results[0,index] = mean_squared_error(all_Y_real,all_Y_pred)
        results_train[:,index] = regressor.predict(Xtrain) 
        results_test[:,index] = regressor.predict(Xtest) 
        r2_train = regressor.score(Xtrain,Ytrain)
        r2_test  = regressor.score(Xtest,Ytest)
        if (r2_train < 0):
            r2_train = np.abs(r2_train)
        if (r2_test < 0):
            r2_test = np.abs(r2_test)
        if (r2_train > 1):
            r2_train = 1/r2_train
        if (r2_test > 1):
            r2_test = 1/r2_test

        results_R2_train[0,index] = r2_train
        results_R2_test[0,index]  = r2_test 
        importances = permutation_importance(regressor,Xtrain,Ytrain)
        importances_results[:,index] = importances["importances_mean"]
        print(str(regressor).split("()")[0]+" complete.")
    print("saving results")
    Save_results(Ytrain,results_train,Ytest,results_test,results_R2_train,results_R2_test,MAE_results,MAPE_results,MSE_results,model_names,(polymer_names,number_of_train_data),importances_results,X_fields,"ET_Optimized_predictions")
    
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=TEST_SIZE, random_state = 42)
    Xtrain = np.array(Xtrain)
    Xtest  = np.array(Xtest)
    Ytrain = np.array(Ytrain).reshape(len(Ytrain),1)
    Ytest  = np.array(Ytest).reshape(len(Ytest),1)
    model_names = [str(name).split("()")[0] for name in GPR_models]
    number_of_features = X.shape[1]
    number_of_models = len(model_names)
    all_Y_real = np.concatenate((Ytrain,Ytest),axis=0) # for calculating MAE, MAPE and MSE
    index = 0
    for Tg in all_Y_real.flatten():
        # find the index of the polymer with the matching Tg value
        Tg_index = polymer_names_with_Tg[polymer_names_with_Tg["nano_composite_Tg_(K)"] == Tg].index
        if not Tg_index.empty:  # check if any index was found
            polymer_names[index] = polymer_names_with_Tg.iloc[Tg_index[0]]["polymer_name"]
            index += 1
        else:
            print("could not find this Tg",Tg)

    len_Xtest  = Xtest.shape[0]
    len_Xtrain = Xtrain.shape[0]
    results_train = np.zeros((len_Xtrain,number_of_models))
    results_test = np.zeros((len_Xtest,number_of_models))
    results_R2_train = np.zeros((1,number_of_models))
    results_R2_test = np.zeros((1,number_of_models))
    importances_results = np.zeros((number_of_features,number_of_models))
    MAE_results = np.zeros((1,number_of_models))
    MAPE_results = np.zeros((1,number_of_models))
    MSE_results = np.zeros((1,number_of_models))
    number_of_train_data = Ytrain.shape[0]
    
    for index, regressor in enumerate(GPR_models):    
        regressor.fit(Xtrain,Ytrain.ravel())
        predict_train = regressor.predict(Xtrain)
        predict_test = regressor.predict(Xtest)
        all_Y_pred = np.concatenate((predict_train,predict_test))
        MAE_results[0,index] = mean_absolute_error(all_Y_real,all_Y_pred)
        MAPE_results[0,index] = mean_absolute_percentage_error(all_Y_real,all_Y_pred)
        MSE_results[0,index] = mean_squared_error(all_Y_real,all_Y_pred)
        results_train[:,index] = regressor.predict(Xtrain) 
        results_test[:,index] = regressor.predict(Xtest) 
        r2_train = regressor.score(Xtrain,Ytrain)
        r2_test  = regressor.score(Xtest,Ytest)
        if (r2_train < 0):
            r2_train = np.abs(r2_train)
        if (r2_test < 0):
            r2_test = np.abs(r2_test)
        if (r2_train > 1):
            r2_train = 1/r2_train
        if (r2_test > 1):
            r2_test = 1/r2_test
        results_R2_train[0,index] = r2_train
        results_R2_test[0,index]  = r2_test 
        importances = permutation_importance(regressor,Xtrain,Ytrain)
        importances_results[:,index] = importances["importances_mean"]
        print(str(regressor).split("()")[0]+" complete.")
    print("saving results")
    Save_results(Ytrain,results_train,Ytest,results_test,results_R2_train,results_R2_test,MAE_results,MAPE_results,MSE_results,model_names,(polymer_names,number_of_train_data),importances_results,X_fields,"GPR_Optimized_predictions")

def Save_results(Ytrain,results_train,Ytest,results_test,results_R2_train,results_R2_test,MAE_results,MAPE_results,MSE_results,model_names,polymer_names_info,importances_result,descriptor_names,filename="results"):
    polymer_names, number_of_train_data = polymer_names_info
    workbook = xlsxwriter.Workbook("results/"+filename+".xlsx")
    worksheet_results = workbook.add_worksheet("results")
    worksheet_metrics = workbook.add_worksheet("metrics")
    worksheet_importances = workbook.add_worksheet("importances")
    
    forward_results         = 0
    forward_R2s             = 0
    title_format = workbook.add_format()
    title_format.set_bold()
    title_format.set_align("center")
    title_format.set_align("vcenter")
    # format for other cells 
    normal_format = workbook.add_format()
    normal_format.set_align("center")
    normal_format.set_align("vcenter")
    # name format
    name_format = workbook.add_format()
    name_format.set_bold()
    name_format.set_align("center")
    name_format.set_align("vcenter")
    name_format.set_pattern(1)
    #name_format.set_bg_color("#89A7DE")
    name_format.set_fg_color("#89A7DE")

    worksheet_importances.write_column(1,0,descriptor_names)
    worksheet_metrics.write(0,0,"names",title_format)
    worksheet_metrics.write(1,0,"R2 train" ,title_format)
    worksheet_metrics.write(2,0,"R2 test" ,title_format)
    worksheet_metrics.write(3,0,"MAE"  ,title_format)
    worksheet_metrics.write(4,0,"MAPE"  ,title_format)
    worksheet_metrics.write(5,0,"MSE"  ,title_format)
    
    for index, name in enumerate(model_names):
    # results:
        worksheet_results.merge_range(0,10*index,0,10*index+9,name,name_format)
        #worksheet_results.write(1,index+forward_results,"R2 train: " + str(results_R2_train[0,index]),normal_format)
        #worksheet_results.write(2,index+forward_results,"R2 test: "  + str(results_R2_test[0,index]),normal_format)
        worksheet_results.write(1,10*index,"polymer names",title_format)
        worksheet_results.write_column(2,10*index,polymer_names[0:number_of_train_data],normal_format)
        worksheet_results.write(1,10*index+1,"Ytrain",title_format)
        worksheet_results.write_column(2,10*index+1,Ytrain,normal_format)
        worksheet_results.write(1,10*index+2,"prediction_train",title_format)
        worksheet_results.write_column(2,10*index+2,results_train[:,index],normal_format)
        worksheet_results.write(1,10*index+3,"MAE",title_format)
        worksheet_results.write_column(2,10*index+3,np.abs(Ytrain.ravel() - results_train[:,index]),normal_format)
        worksheet_results.write(1,10*index+4,"MAPE",title_format)
        worksheet_results.write_column(2,10*index+4,np.divide(np.abs(Ytrain.ravel() - results_train[:,index])*100,Ytrain.ravel()),normal_format)

        worksheet_results.write(1,10*index+5,"polymer names",title_format)
        worksheet_results.write_column(2,10*index+5,polymer_names[number_of_train_data:],normal_format)
        worksheet_results.write(1,10*index+6,"Ytest",title_format)
        worksheet_results.write_column(2,10*index+6,Ytest,normal_format)
        worksheet_results.write(1,10*index+7,"prediction_test",title_format)
        worksheet_results.write_column(2,10*index+7,results_test[:,index],normal_format)
        worksheet_results.write(1,10*index+8,"MAE",title_format)
        worksheet_results.write_column(2,10*index+8,np.abs(Ytest.ravel() - results_test[:,index]),normal_format)
        worksheet_results.write(1,10*index+9,"MAPE",title_format)
        worksheet_results.write_column(2,10*index+9,np.divide(np.abs(Ytest.ravel() - results_test[:,index])*100,Ytest.ravel()),normal_format)
    # metrics:
        worksheet_metrics.write(0,index+1,name,title_format)
        worksheet_metrics.write(1,index+1,results_R2_train[0,index],normal_format)
        worksheet_metrics.write(2,index+1,results_R2_test[0,index] ,normal_format)
        worksheet_metrics.write(3,index+1,MAE_results[0,index] ,normal_format)
        worksheet_metrics.write(4,index+1,MAPE_results[0,index] *100,normal_format)
        worksheet_metrics.write(5,index+1,MSE_results[0,index] ,normal_format)

    # importances
        worksheet_importances.write(0,0, "name",title_format)
        worksheet_importances.write(0,index+1,name,title_format)
        worksheet_importances.write_column(1,index+1,importances_result[1:,index],normal_format)
        forward_results += 12
        forward_R2s += 1
    worksheet_results.autofit()
    worksheet_metrics.autofit()
    worksheet_importances.autofit()
    workbook.close()

results = []

Optimized_payan_name_Tg_nano_composite()