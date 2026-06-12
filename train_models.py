import pandas as pd
import xlsxwriter
from sklearn.ensemble import GradientBoostingRegressor ,ExtraTreesRegressor,RandomForestRegressor 
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, Matern, Exponentiation
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error, mean_squared_error
import numpy as np
import pickle as pkl

def Optimized_Tg_nano_composite():
    models = [
            # Extra Trees optimized results
            ExtraTreesRegressor(n_estimators = 272, min_samples_split = 0.089794412),
            ExtraTreesRegressor(n_estimators = 172, min_samples_split = 0.134436839),
            ExtraTreesRegressor(n_estimators = 224, min_samples_split = 0.143697764),
            ExtraTreesRegressor(n_estimators = 153, min_samples_split = 0.075041656),
            ExtraTreesRegressor(n_estimators = 172, min_samples_split = 0.134436839),
            # Random Forest
            RandomForestRegressor(n_estimators = 335 , min_samples_split = 0.06259375 ),
            RandomForestRegressor(n_estimators = 334 , min_samples_split = 0.089794412),
            RandomForestRegressor(n_estimators = 100 , min_samples_split = 0.065090554),
            RandomForestRegressor(n_estimators = 311 , min_samples_split = 0.057541194),
            RandomForestRegressor(n_estimators = 334 , min_samples_split = 0.089794412),
            # Gradient Boosting optimized results
            GradientBoostingRegressor(n_estimators = 178, min_samples_split = 0.693385724),
            GradientBoostingRegressor(n_estimators = 193, min_samples_split = 0.705201292),
            GradientBoostingRegressor(n_estimators = 193, min_samples_split = 0.705201292),
            GradientBoostingRegressor(n_estimators = 186, min_samples_split = 0.69798037 ),
            GradientBoostingRegressor(n_estimators = 193, min_samples_split = 0.705201292),
            # Gaussian Process RBF optimized results
            GaussianProcessRegressor(kernel = RBF(length_scale = 1.617183123), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel = RBF(length_scale = 8.45381705 ), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel = RBF(length_scale = 2.939111441), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel = RBF(length_scale = 0.012734077), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel = RBF(length_scale = 8.45381705 ), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            #Gaussian Process Exponent RBF optimized results
            GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 0.504308474),exponent = 0.093200173), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 3.926932014),exponent = 0.051818176), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 3.373040455),exponent = 0.070473753), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 5.827752653),exponent = 2.979309401), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel=Exponentiation(kernel=RBF(length_scale = 3.926932014),exponent = 0.051818176), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            # Gaussian Process Matern optimized results
            GaussianProcessRegressor(kernel=Matern(length_scale = 1.127365022, nu = 0.5), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel=Matern(length_scale = 0.192552148, nu = 0.5), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel=Matern(length_scale = 4.570861238, nu = 0.5), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel=Matern(length_scale = 4.638781749, nu = 0.5), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel=Matern(length_scale = 0.192552148, nu = 0.5), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            # Gaussian Process Exponent Matern optimized results
            GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 1.545304196, nu = 0.5), exponent = 1.346181932), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 0.592311091, nu = 0.5), exponent = 1.470672505), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 0.247757204, nu = 0.5), exponent = 1.033802531), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 1.132346208, nu = 0.5), exponent = 1.428879482), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5),
            GaussianProcessRegressor(kernel = Exponentiation(kernel = Matern(length_scale = 0.592311091, nu = 0.5), exponent = 1.470672505), normalize_y = True, alpha=1e-4, random_state = 42,  n_restarts_optimizer=5)
        ]
    
    filename = "dataset"
    dataset = pkl.load(open("data/"+filename+".rb", "rb"))

    X_fields = ["graphene_w","CH","C=OH","OH","NH","B","BCH","hydrogen_bond","polarity"]
    Y_field = ["nano_composite_Tg_(K)"]
    X = dataset.iloc[:,[4,9,11,15,17,19,21,35,36,37]]
    Y = dataset.iloc[:,5]
    TEST_SIZE = 0.15
    polymer_names_with_Tg = dataset.iloc[:,[0,5]]
    
    # training models
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=TEST_SIZE, random_state = 42)
    scaler = StandardScaler()
    Xtrain = scaler.fit_transform(Xtrain)
    Xtest  = scaler.transform(Xtest)
    Ytrain = np.array(Ytrain).reshape(len(Ytrain),1)
    Ytest  = np.array(Ytest).reshape(len(Ytest),1)
    model_names = [str(name).split("()")[0] for name in models]
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
    
    for index, regressor in enumerate(models):    
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
    Save_results(Ytrain,results_train,Ytest,results_test,results_R2_train,results_R2_test,MAE_results,MAPE_results,MSE_results,model_names,(polymer_names,number_of_train_data),importances_results,X_fields,"Optimized_predictions")
    
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

Optimized_Tg_nano_composite()