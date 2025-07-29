classdef MethodsFigureTimeWindow < MethodsFigure

    % created in order to clone methodsFigure for a version with a specific
    % timewindow.

    % See MethodsFigure and FigurePipeline for included functions. 



    properties
    end


    methods

        function obj = MethodsFigureTimeWindow(MethodsFigure,timeWindow)
            %timeWindow must be created within MethodsFigure.
            %To see implementation of methodsFigure, see FigureScript (RW)

            obj.positionData = MethodsFigure.positionData.getTimeWindow(timeWindow);
           % obj.positionDataLinearized = MethodsFigure.positionDataLinearized.getTimeWindow(timeWindow);
            obj.LFP = MethodsFigure.LFP.getTimeWindow(timeWindow);

        end




    end













end