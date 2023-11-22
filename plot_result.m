% Copyright 2023, Chuan Huang
    
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at

%   http://www.apache.org/licenses/LICENSE-2.0

% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
function [h] = plot_result(output)

    %% plot
    h=tiledlayout(4, 2);
    nexttile([2, 2]);
    plot(output.GT.pos(1,:), output.GT.pos(2,:),'m', 'LineWidth', 2)
    hold on;
    plot(output.INS.pos(1,:), output.INS.pos(2,:),'b', 'LineWidth', 2)
    plot(output.MAINS.pos(1,:), output.MAINS.pos(2,:),'k', 'LineWidth', 2)
    plot(output.hassen.pos(1,:), output.hassen.pos(2,:),'Color',[0.8500 0.3250 0.0980], 'LineWidth', 2)

    t = (10: length(output.INS.pos) - 1) * 0.01;
    legend('True trajectory', 'INS trajectory', 'MAINS trajectory' , 'Trajectory from [19]','FontName','Times New Roman','FontSize',14)
    title('Horizontal trajectory','FontName','Times New Roman','FontSize',18)
    xlabel('x [m]','FontName','Times New Roman','FontSize',17)
    ylabel('y [m]','FontName','Times New Roman','FontSize',17)
    ax = gca;
    ax.XAxis.FontSize = 14; 
    ax.YAxis.FontSize = 14;
    xlim([-12 12])
    ylim([-12 12])
    grid on;
    nexttile([1, 2]);
    semilogy(t, output.MAINS.err(1,11:end),'k', 'LineWidth', 2)
    hold on
    semilogy(t, output.INS.err(1,11:end),'b', 'LineWidth', 2)
    semilogy(t, output.hassen.err(1,11:end),'Color',[0.8500 0.3250 0.0980], 'LineWidth', 2)
    title('Horizontal error','FontName','Times New Roman','FontSize',18)
    xlabel('time [s]','FontName','Times New Roman','FontSize',17)
    ylabel('error [m]','FontName','Times New Roman','FontSize',17)
    ylim([0, 2e2]);
    yticks([1e-5 1 ])
    ax = gca;
    ax.XAxis.FontSize = 14; 
    ax.YAxis.FontSize = 14;
    grid on;
    nexttile([1, 2]);
    semilogy(t, output.MAINS.err(2,11:end),'k', 'LineWidth', 2)
    hold on
    semilogy(t, output.INS.err(2,11:end),'b', 'LineWidth', 2)
    semilogy(t, output.hassen.err(2,11:end),'Color',[0.8500 0.3250 0.0980], 'LineWidth', 2)
    title('Vertical error','FontName','Times New Roman','FontSize',18)
    xlabel('time [s]','FontName','Times New Roman','FontSize',17)
    ylabel('error [m]','FontName','Times New Roman','FontSize',17)
    ylim([0, 1e2]);
    yticks([1e-5 1 ])
    ax = gca;
    ax.XAxis.FontSize = 14; 
    ax.YAxis.FontSize = 14;
    grid on;
    h.Padding='compact';


end