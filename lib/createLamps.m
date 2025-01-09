function [fig] = createLamps(x,chann_names)
    
    ScrnSz= get(0,"ScreenSize");
    ScreenHeight = ScrnSz(4);%*system_scaling;
    disp(ScreenHeight)
    ScreenWidth = ScrnSz(3);%*system_scaling;
    disp(ScreenWidth)

    % Define lamp size, spacing, and initial position
    lampSize = [16, 16];  % width, height
    spacing = 24;
    initPosition = [10, ScreenHeight-50];  % starting position for the first lamp
    max_lamp_1col = floor((ScreenHeight-30-50)/spacing);
    % Create lamps
    label_width = 40;
    msg_width = 25;
    

    
    intercol_gap = label_width + lampSize(1) + spacing + msg_width;


    %create a figure
    wd_width = 120 * ceil((x+4)/max_lamp_1col);

    fig = figure('Name', 'Quality Indicators', 'NumberTitle', 'off','ToolBar', 'none','MenuBar','none', ...
                 'Position', [ScreenWidth-wd_width, 0, wd_width, ScreenHeight]); % [left, bottom, width, height]
    
    fig.UserData.lamp_handles = gobjects(1, x); % Preallocate array for lamp handles
    fig.UserData.msg_handles = gobjects(1, x); % Preallocate array for text handles
    
    delete(fig.Children);
    colTitles = {"CHANN", "", "CAUSE"};

    %set(fig, 'WindowStyle', 'alwaysontop'); % Make the figure always on top

    col_pos = {[initPosition(1),initPosition(2),label_width,lampSize(2)],[initPosition(1)+label_width,initPosition(2),label_width,lampSize(2)],...
        [initPosition(1)+label_width+msg_width,initPosition(2),label_width,lampSize(2)]};
    for col = 1:3
        uicontrol(fig, 'Style', 'text', 'Position',col_pos{col}, ...
            'String', colTitles{col}, 'FontSize', 6, 'FontWeight', 'bold', ...
            'BackgroundColor', fig.Color, 'HorizontalAlignment', 'left');
    end
    for i = 1:x
        

        curr_col = ceil(i/max_lamp_1col);
        curr_row = mod(i-1, max_lamp_1col);
        l_pos = initPosition(1) + (curr_col-1)*intercol_gap; % left position
        b_pos = initPosition(2) - (curr_row+1) * spacing; % bottom position

        labelPosition = [l_pos, b_pos, label_width, lampSize(2)];% [left, bottom, width, height];
        lamp_position = [l_pos + label_width, b_pos, lampSize(1), lampSize(2)];% [left, bottom, width, height]
        msg_position = [l_pos + label_width+msg_width, b_pos, lampSize(1)+12, lampSize(2)];% [left, bottom, width, height]
        fig.UserData.lamp_handles(i) = uicontrol('Style', 'pushbutton', 'Position', lamp_position, 'BackgroundColor', 'green', 'Enable', 'inactive');
        uicontrol(fig, 'Style', 'text', 'Position', labelPosition, 'String', chann_names{i},'FontSize', 10, 'BackgroundColor', fig.Color, 'HorizontalAlignment', 'left');
        fig.UserData.msg_handles(i) =  uicontrol(fig, 'Style', 'text', ...
                                            'Position', msg_position, ...
                                            'String', "", ...
                                            'FontSize', 10, ...
                                            'BackgroundColor', fig.Color, ...
                                            'HorizontalAlignment', 'left');

        % if i == x && b_pos > 6*spacing
        %     for j = 1:5
        %         uicontrol(fig, 'Style', 'text', 'Position',[l_pos,5+j*spacing,intercol_gap,lampSize(2)], ...
        %           'String', 'MU: Multiple Causes', 'FontSize', 10, ...
        %           'BackgroundColor', 'blue', 'HorizontalAlignment', 'left');
        %     end
        % end
    end

    if i == x && b_pos > 6*spacing
    % Define panel size and position
        panelWidth = intercol_gap+10;
        panelHeight = 5 * spacing + lampSize(2);
        panelLeft = l_pos;
        panelBottom = 35; % Bottom position for the panel
       
    else
        panelWidth = intercol_gap+10;
        panelHeight = 5 * spacing + lampSize(2);
        panelLeft = l_pos+intercol_gap;
        panelBottom = 35; % Bottom position for the panel
    end
    % Create a panel with a background color
    legendPanel = uipanel(fig, 'Position', [panelLeft/fig.Position(3), panelBottom/fig.Position(4), ...
        panelWidth/fig.Position(3), panelHeight/fig.Position(4)], ...
        'BackgroundColor', 'blue', 'BorderType', 'etchedin');


    % Add legend entries inside the panel
    legendEntries = {'MU: Multiple Causes', 'NS: No Signal', 'HA: High Amplitude', ...
        'FN: Frequency Noise', 'LC: Low Correlation'};
    for j = 1:numel(legendEntries)
        uicontrol(legendPanel, 'Style', 'text', 'Position', [5, panelHeight - j*spacing, panelWidth+10, lampSize(2)], ...
            'String', legendEntries{j}, 'FontSize', 8,...
            'BackgroundColor', 'blue', 'ForegroundColor', 'white', 'HorizontalAlignment', 'left');
    end

    

end

