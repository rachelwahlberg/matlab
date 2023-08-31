classdef RWexperiment < Singleton
    %RWEXPERIMENT - Rachel Wahlberg adaptation of utku's SDExperiment class.
    %   Detailed explanation goes here
    
    properties
        XMLFile
    end
    methods(Access=private)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function newObj = RWexperiment(xmlfile)
            % Initialise your custom properties.
            if exist("xmlfile",'var')
                newObj.XMLFile = xmlfile;
            else
            newObj.XMLFile = ...
                './ExperimentSpecific/Configure/Experiment.xml';
            end
            try
                S=readstruct(newObj.XMLFile);
            catch
                S.FileLocations.General.ExperimentConfig=newObj.XMLFile;
                writestruct(S,newObj.XMLFile)
            end
            try structstruct(S); catch, end
        end
    end
    
    methods(Static)
        % Concrete implementation.  See Singleton superclass.
        function obj = instance(xmlFile)
            persistent uniqueInstance
            if isempty(uniqueInstance)
                if exist('xmlFile','var')
                    obj = experiment.SDExperiment(xmlFile);
                else
                    obj = experiment.SDExperiment();
                end
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    end
    
    
    methods
        function S = get(obj)
            S=readstruct(obj.XMLFile);
        end
        function color = getStateColors(obj,state)
            S=readstruct(obj.XMLFile);
            statestr=fieldnames(S.Colors);
            for ist=1:numel(statestr)
                statecodes(ist)=S.StateCodes.(statestr{ist});
            end
            for ist=1:numel(statestr)
                colors{ist}=S.Colors.(statestr{ist})/255;
            end
            color=containers.Map(statecodes, colors);
            if exist('state','var')
                if isnumeric( state)
                    statecode=state;
                else
                    statecode=S.StateCodes.(state);
                end
                color=color(statecode);
            end
        end
        function state1 = getStateCode(obj,state)
            S=readstruct(obj.XMLFile);
            sc=fieldnames(   S.StateCodes);
            for ist=1:numel(sc)
                statecodes(ist)=S.StateCodes.(sc{ist}); %#ok<AGROW>
            end
            states=containers.Map(statecodes, sc);
           
            if exist('state','var')
                if isnumeric( state)
                    state1=states(state);
                    
                else
                    idx=ismember(sc,state);
                    state1=statecodes(idx);
                end
            end
        end
        function S = set(obj,S)
            writestruct(S,obj.XMLFile);
            S=obj.get;
            try structstruct(S); catch, end
        end
    end
end

