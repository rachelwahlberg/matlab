function parentname = initializeDirectories(basepath,rat,sessiontype,region,folderdate)

if ~isfolder([basepath '\' rat])
   mkdir([basepath '\' rat]);
end

if strcmp(sessiontype,'Recording') == 1

    parentname = [basepath '\' rat '\Recordings\' region '\' folderdate];

    if ~isfolder(parentname) %if ~isfolder([rats_basepath '\' rat '\**\' currentdate])

    mkdir(parentname)
    mkdir([parentname '\OpenEphys'])
    mkdir([parentname '\National_Instruments'])
    mkdir([parentname '\Videos'])
    mkdir([parentname '\Figures'])
    mkdir([parentname '\Motive'])
    end

elseif strcmp(sessiontype,'TaskTraining') == 1

    parentname = [basepath '\' rat '\TaskTraining\' folderdate];

    if ~isfolder(parentname)

    mkdir(parentname)
    mkdir([parentname '\National_Instruments'])
    mkdir([parentname '\Videos'])
    mkdir([parentname '\Figures'])

    end

elseif strcmp(sessiontype,'WaterTraining') == 1

    parentname = [basepath '\' rat '\WaterTraining\' folderdate];

    if ~isfolder(parentname)

    mkdir(parentname)
    mkdir([parentname '\National_Instruments'])
    mkdir([parentname '\Videos'])
    mkdir([parentname '\Figures'])

    end

elseif strcmp(sessiontype,'TEST') == 1

    parentname = [basepath '\' rat '\TEST\' folderdate];

    if ~isfolder(parentname)

    mkdir(parentname)
    mkdir([parentname '\National_Instruments'])
    mkdir([parentname '\Videos'])
    mkdir([parentname '\Figures'])

    end

else
   disp('oops check your sessiontype name') 
   return
end







































end