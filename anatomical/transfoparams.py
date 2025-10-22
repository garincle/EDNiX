
default_intramod_type = 'SyNCC' # should be "SyN","SyNCC","SynBold","Rigid" or "Affine"
default_intermod_type = 'SyN'   # should be "SyN","SyNCC","SynBold","Rigid" or "Affine"
default_affmetric     = 'MI'    # should be "mattes" (means MI for mutual information),"GC" (global correlation) or "meansquare"
default_affmetricT    = 'mattes' # should be "mattes" (means MI for mutual information),"GC" (global correlation) or "meansquare"
default_interp        = 'hammingWindowedSinc' # should be 'nearestNeighbour','linear','hammingWindowedSinc'


def set(transfo_message,IgotbothT1T2,creat_study_template,brain_skullstrip_1,brain_skullstrip_2,template_skullstrip,Align_img_to_template,list_transfo):

    if transfo_message == 'do_it_for_me':
        list_transfo = []
        # presume that every parameter is empty

        Align_img_to_template = 'Ants'  # should be "3dAllineate" / "No" / "@Align_Centers" / " Ants"

        transfo_dict = {"name": 'align',
                              "interpol": default_interp,
                              "type_of_transform": 'Rigid',
                              "affmetric": default_affmetric,
                              "affmetricT": default_affmetricT
                              }
        list_transfo.append(transfo_dict)

        if IgotbothT1T2 == True:
            # specific for intermodality
            transfo_dict["name"] = 'between'
            list_transfo.append(transfo_dict)

        transfo_dict = {"name": 'coreg',"interpol": default_interp,
                                "type_of_transform": default_intramod_type,
                                "affmetric": default_affmetric,
                                "affmetricT": default_affmetricT
                                }
        list_transfo.append(transfo_dict)

        check_ants1 = brain_skullstrip_1.split('_')
        if 'ANTS' in check_ants1:
            transfo_dict["name"] = 'SS1'
            list_transfo.append(transfo_dict)

        check_ants2 = brain_skullstrip_2.split('_')
        if 'ANTS' in check_ants2:
            transfo_dict["name"] = 'SS2'
            list_transfo.append(transfo_dict)

        if creat_study_template == True:
            check_ants3 = template_skullstrip.split('_')
            if 'ANTS' in check_ants3:
                transfo_dict["name"] = 'SS3'
                list_transfo.append(transfo_dict)
            transfo_dict["name"] = 'stdyT'
            list_transfo.append(transfo_dict)



    if transfo_message == 'do_as_I_said':

        if Align_img_to_template == '':
            Align_img_to_template = 'Ants'
            for i, j in enumerate(list_transfo):
                if list_transfo[i]["name"] == 'align':
                    list_transfo[i]["interpol"]          = default_interp
                    list_transfo[i]["type_of_transform"] = 'Rigid'
                    list_transfo[i]["affmetric"]         = default_affmetric
                    list_transfo[i]["affmetricT"]        = default_affmetricT

        refnb = 0
        for i, j in enumerate(list_transfo):
            if list_transfo[i]["name"] == 'coreg':
                refnb = i

        dummy_list = ['SS1','SS2','SS3','stdyT']


        for i, j in enumerate(list_transfo):
            if list_transfo[i]["name"] == 'align' or list_transfo[i]["name"] == 'between':
                for default,dictname in zip([default_interp,'Rigid',default_affmetric,default_affmetricT],
                                        ["interpol","type_of_transform","affmetric","affmetricT"]):
                    if list_transfo[i][dictname] == '' :
                        list_transfo[i][dictname] = default

            if list_transfo[i]["name"] == 'coreg':
                for default,dictname in zip([default_interp,default_intramod_type,default_affmetric,default_affmetricT],
                                        ["interpol","type_of_transform","affmetric","affmetricT"]):
                    if list_transfo[i][dictname] == '' :
                        list_transfo[i][dictname] = default

            if any(list_transfo[i]["name"] in word for word in dummy_list):
                for default,dictname in zip([default_interp,default_intramod_type,default_affmetric,default_affmetricT],
                                        ["interpol","type_of_transform","affmetric","affmetricT"]):
                    if list_transfo[i][dictname] == '' :
                        list_transfo[i][dictname] = list_transfo[refnb][dictname]












    return (Align_img_to_template,list_transfo)


