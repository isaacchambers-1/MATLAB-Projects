classdef WheelNoiseApp < handle
    % WheelNoiseApp – Configure wheel-speed sensor noise over a UDDS drive cycle.
    % Generates a *masked*, unevenly-spaced data set when dropout occurs
    % (no filling). Exports/assigns masked vectors for student analysis.
    %
    % Output table S:
    %   time  : masked time vector (dropout samples removed)
    %   v_true: true UDDS speed at those same times (masked)
    %   v_wheel: noisy wheel-speed samples (masked)
    %
    % Buttons:
    %   • Generate
    %   • Send S to Workspace
    %   • Export CSVs (UDDS_true.csv, UDDS_wheel.csv)  <-- masked, same length

    properties (Access = private)
        UIF          matlab.ui.Figure
        Grid         matlab.ui.container.GridLayout
        Left         matlab.ui.container.Panel
        Right        matlab.ui.container.Panel
        GridCtl      matlab.ui.container.GridLayout
        % Controls
        fileField    matlab.ui.control.EditField
        browseBtn    matlab.ui.control.Button
        seedField    matlab.ui.control.NumericEditField
        sigmaField   matlab.ui.control.NumericEditField
        quantField   matlab.ui.control.NumericEditField
        dropField    matlab.ui.control.NumericEditField
        zoomStart    matlab.ui.control.NumericEditField
        zoomEnd      matlab.ui.control.NumericEditField
        genBtn       matlab.ui.control.Button
        toWSBtn      matlab.ui.control.Button
        exportBtn    matlab.ui.control.Button
        % Axes
        AxFull       matlab.ui.control.UIAxes
        AxZoom       matlab.ui.control.UIAxes
        % Data
        T            table
        S            table
        t_full       double   % full (unmasked) time for plotting context
        v_true_full  double   % full (unmasked) UDDS true for plotting context
    end

    methods
        function app = WheelNoiseApp()
            buildUI(app);
            % Try to auto-load UDDS.csv in current folder
            defaultFile = fullfile(pwd,'UDDS.csv');
            if exist(defaultFile,'file')
                app.fileField.Value = defaultFile;
                try
                    app.T = readtable(defaultFile);
                catch
                    app.T = table();
                end
            end
        end
    end

    methods (Access = private)
        function buildUI(app)
            app.UIF = uifigure('Name','Wheel-Speed Noise Generator','Position',[100 100 1100 680]);
            app.Grid = uigridlayout(app.UIF,[1 2]); app.Grid.ColumnWidth = {360,'1x'};

            % Left panel (controls)
            app.Left = uipanel(app.Grid,'Title','Parameters');
            app.GridCtl = uigridlayout(app.Left,[12 2]);
            app.GridCtl.RowHeight = repmat({'fit'},1,12);
            app.GridCtl.ColumnWidth = {'1x','1x'};

            % File & seed
            uilabel(app.GridCtl,'Text','UDDS file');
            fileGrid = uigridlayout(app.GridCtl,[1 2]); fileGrid.ColumnWidth={'1x',70};
            app.fileField = uieditfield(fileGrid,'text','Value','UDDS.csv');
            app.browseBtn = uibutton(fileGrid,'Text','Browse','ButtonPushedFcn',@(s,e)onBrowse(app));

            uilabel(app.GridCtl,'Text','Seed');
            app.seedField = uieditfield(app.GridCtl,'numeric','Value',42,'Limits',[-Inf Inf]);

            % Noise controls
            uilabel(app.GridCtl,'Text','Gaussian σ [m/s]');
            app.sigmaField = uieditfield(app.GridCtl,'numeric','Value',0.02,'Limits',[0 Inf]);

            uilabel(app.GridCtl,'Text','Quantization [m/s]');
            app.quantField = uieditfield(app.GridCtl,'numeric','Value',0.05,'Limits',[0 Inf]);

            uilabel(app.GridCtl,'Text','Dropout frac [0..1]');
            app.dropField  = uieditfield(app.GridCtl,'numeric','Value',0.005,'Limits',[0 1]);

            % Zoom window
            uilabel(app.GridCtl,'Text','Zoom start [s]');
            app.zoomStart = uieditfield(app.GridCtl,'numeric','Value',0,'Limits',[-Inf Inf]);
            uilabel(app.GridCtl,'Text','Zoom end [s]');
            app.zoomEnd   = uieditfield(app.GridCtl,'numeric','Value',0,'Limits',[-Inf Inf]);

            % Buttons
            app.genBtn   = uibutton(app.GridCtl,'Text','Generate','ButtonPushedFcn',@(s,e)onGenerate(app));
            app.toWSBtn  = uibutton(app.GridCtl,'Text','Send S to Workspace','ButtonPushedFcn',@(s,e)onToWS(app));
            app.exportBtn= uibutton(app.GridCtl,'Text','Export CSVs','ButtonPushedFcn',@(s,e)onExport(app));

            % Right panel (plots)
            app.Right = uipanel(app.Grid,'Title','Plots');
            rightGrid = uigridlayout(app.Right,[2 1]); rightGrid.RowHeight = {'1x','1x'};
            app.AxFull = uiaxes(rightGrid); title(app.AxFull,'Full UDDS (context)'); ylabel(app.AxFull,'Speed [m/s]'); grid(app.AxFull,'on');
            app.AxZoom = uiaxes(rightGrid); title(app.AxZoom,'Zoomed Segment'); xlabel(app.AxZoom,'Time [s]'); ylabel(app.AxZoom,'Speed [m/s]'); grid(app.AxZoom,'on');
        end

        function onBrowse(app)
            [f,p] = uigetfile({'*.csv;*.txt','Data files (*.csv,*.txt)'},'Select UDDS file');
            if isequal(f,0), return; end
            app.fileField.Value = fullfile(p,f);
            try
                app.T = readtable(fullfile(p,f));
                uialert(app.UIF,'File loaded. Click Generate to synthesize wheel-speed signal.','Loaded');
            catch ME
                uialert(app.UIF,ME.message,'Load Error');
            end
        end

        function onGenerate(app)
            % Load/refresh table
            try
                app.T = readtable(app.fileField.Value);
            catch ME
                uialert(app.UIF,['Failed to read file: ' ME.message],'Read error'); return;
            end

            % Extract time/speed (full, unmasked) + store for context plots
            [t_full, v_true_full] = WheelNoiseApp.load_udds_like(app.T);
            app.t_full = t_full; app.v_true_full = v_true_full;

            % RNG seed
            if ~isfinite(app.seedField.Value), rng('default'); else, rng(app.seedField.Value); end

            % Noise params
            sigma = app.sigmaField.Value;
            q     = app.quantField.Value;
            pdrop = app.dropField.Value;

            % Noisy wheel (introduces NaNs for dropout)
            v_wheel_full = WheelNoiseApp.add_noise_core(v_true_full, sigma, q, pdrop);

            % === MASK & CREATE UNEVENLY-SPACED VECTORS (no fill) ===
            keep = ~isnan(v_wheel_full);
            t    = t_full(keep);
            v_w  = v_wheel_full(keep);
            v_tr = v_true_full(keep);  % true at the same surviving time stamps

            % Package masked table
            app.S = table(t, v_tr, v_w, 'VariableNames',{'time','v_true','v_wheel'});

            % -------- Plots --------
            % Full: show full true (gray), masked true + masked wheel
            cla(app.AxFull);
            plot(app.AxFull, t_full, v_true_full, 'DisplayName','UDDS true (full)'); hold(app.AxFull,'on');
            plot(app.AxFull, t, v_tr,  '-',  'DisplayName','UDDS true (masked)');
            plot(app.AxFull, t, v_w,   '--',   'DisplayName','Wheel (noisy, masked)');
            legend(app.AxFull,'Location','best'); hold(app.AxFull,'off');

            % Zoom window relative to full timebase
            if isfinite(app.zoomStart.Value) && isfinite(app.zoomEnd.Value) && ...
               (app.zoomEnd.Value > app.zoomStart.Value) && ...
               ~(app.zoomStart.Value==0 && app.zoomEnd.Value==0)
                t1 = app.zoomStart.Value; t2 = app.zoomEnd.Value;
            else
                mid = 0.5*(t_full(1)+t_full(end)); t1 = max(t_full(1),mid-30); t2 = min(t_full(end),mid+30);
            end
            idx_full = t_full>=t1 & app.t_full<=t2;
            idx_mask = t>=t1 & t<=t2;

            cla(app.AxZoom);
            plot(app.AxZoom, t_full(idx_full), v_true_full(idx_full), 'DisplayName','UDDS true (full)'); hold(app.AxZoom,'on');
            plot(app.AxZoom, t(idx_mask), v_tr(idx_mask), '-', 'DisplayName','UDDS true (masked)');
            plot(app.AxZoom, t(idx_mask), v_w(idx_mask),  '--',  'DisplayName','Wheel (noisy, masked)');
            title(app.AxZoom, sprintf('Zoom: [%.0f, %.0f] s', t1, t2));
            legend(app.AxZoom,'Location','best'); hold(app.AxZoom,'off');

            uialert(app.UIF,'Wheel-speed signal generated (masked, uneven spacing).','Done');
        end

        function onToWS(app)
            if isempty(app.S)
                uialert(app.UIF,'No data yet. Click Generate first.','Workspace'); return;
            end
            assignin('base','S',app.S);
            uialert(app.UIF,'Assigned masked table S to workspace.','Workspace');
        end

        function onExport(app)
            if isempty(app.S)
                uialert(app.UIF,'No data yet. Click Generate first.','Export'); return;
            end
            p = uigetdir(pwd,'Choose export folder');
            if isequal(p,0), return; end
            % Export *masked* vectors so x,y are aligned and unevenly spaced
            Ttrue = table(app.S.time, app.S.v_true,'VariableNames',{'time','v_true'});
            Tw    = table(app.S.time, app.S.v_wheel,'VariableNames',{'time','v_wheel'});
            writetable(Ttrue, fullfile(p,'UDDS_true.csv'));
            writetable(Tw,    fullfile(p,'UDDS_wheel.csv'));
            uialert(app.UIF,'Exported masked UDDS_true.csv and UDDS_wheel.csv.','Export');
        end
    end

    % ------------------------ Static helpers ------------------------
    methods (Access = private, Static)
        function [t,v_true] = load_udds_like(T)
            low = lower(T.Properties.VariableNames);
            tcol = find(ismember(low, {'time','t','time_s'}), 1);
            vcol = find(ismember(low, {'speed','v','speed_mph','velocity'}), 1);
            if isempty(tcol) || isempty(vcol)
                error('Table must contain recognizable time and speed columns.');
            end
            t = T{:,tcol}; t = t(:);
            v = T{:,vcol}; v = v(:);
            if t(1)~=0, t = t - t(1); end
            mu = mean(v,'omitnan');
            if mu > 60
                v_true = v/3.6;       % km/h → m/s
            elseif mu > 30
                v_true = v*0.44704;   % mph → m/s
            else
                v_true = v;           % m/s
            end
        end

        function y = add_noise_core(x, sigma, q, pdrop)
            y = x(:);
            if sigma>0, y = y + sigma*randn(size(y)); end
            if q>0,     y = q*round(y/q);             end
            if pdrop>0
                mask = rand(size(y)) < pdrop;
                y(mask) = NaN;
            end
        end
    end
end
