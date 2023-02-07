clc; clear
close all;

%% input arguments to select simulation setup
counter = 0;
while 1
    
    if counter > 0
        prompt1 = 'Continue or not? (y/n): ';
        selection = input(prompt1, 's');
        if selection ~= 'y'
            %             close all;
            %             clc; clear;
            break;
        end
    end
    
    prompt2 = 'Please select a simulation (sim/less/switch/fault): ';
    option = input(prompt2, 's');
    
    switch option
        case {'sim', 'less', 'switch', 'fault'}
            counter = 1;
            
            prompt3 = 'Please select the dynamics (integrator/nonholonomic): ';
            dynamics = input(prompt3, 's');
            
            switch dynamics
                case {'integrator', 'nonholonomic'}
                    % run simulation
                    Ns = input('Please set the number of harmonics: ');
                    N = input('Please set the number of agents: ');
                    [data, sim] = contour(option, dynamics, Ns, N);
                    
                    % gif
                    plotGif (data, Ns, option);
                    % figures
                    plotFigure(data, sim, option);
                    
                otherwise
                    fprintf('%s\n','Wrong arguments! Please try again.' )
            end
            
        otherwise
            counter = 1;
            fprintf('%s\n','Wrong arguments! Please try again.' )
    end
    
end
