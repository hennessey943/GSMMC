function Cell_Growth(A, n, time)
% A is n-by-n matrix
delta = 3; % Scaling Factor
axis( [-delta (n+1)*delta -delta (n+1)*delta] );
for r = 1:n
    for c = 1:n
        s = sqrt( A(r,c) );  % Side Lengths
        x = r*delta; y = c*delta; % Center of Squares
        %if time == 0.5
        %    fill( [x-s/2 x+s/2 x+s/2 x-s/2], [y-s/2 y-s/2 y+s/2 y+s/2], 'b' );
        %    hold on;
        %else
            fill( [x-s/2 x+s/2 x+s/2 x-s/2], [y-s/2 y-s/2 y+s/2 y+s/2], 'g' );
            plot_title = sprintf( '%0.1f days', time );
            title( plot_title );
            hold on;
        %end
        axis( [0 (n+1)*delta 0 (n+1)*delta] );
    end
end
%if time == 0.5
%    pause(1.5);
%end
pause(0.005);
hold off

