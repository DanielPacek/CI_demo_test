clear all; clc; close all;

% DEFINOVANIE DELIACEJ PRIAMKY %   
line_points_1 = [ 0 , 0 ; 0 , 4 ]; % [ Ax , Ay ; Bx , By] 
line_points_2 = [ 8 , 0 ; 11 , 4]; % pre definovanie priamky su pouzite 2 body
plot_line = 9; % pre vykreslenie priamok zadajte 1
                % hodnota ina ako 0 alebo 1 nevykresli priamky
                
global data; % globalna premenna 

data = readmatrix('StepAx2_2019_2_13_ind_mag2.txt'); % nacitanie dat
[length_ , ~] = size(data);
[length_,data] = zero_line_del(data,length_); % vymaze riadky obsahujuce iba nuly

add_reference_number(); %prida referencne cislo,ktore hovori o pocte
                        %udajov pre jedeu polohu Y

% IDENTIFIKACIA TVARU POZADOVANEJ KRIVKY %
ref_2 = 0; 
ref_4 = 0;
ref_6 = 0;
line_1_active = 0; % overi ci bola vyuzita deliaca priamka 1
line_2_active = 0; % overi ci bola vyuzita deliaca priamka 2

for i = 1 : 1 : length_ % zisti ake hodnoty nadobuda ref. cislo
    if data( i , 5 ) == 2
        ref_2 = ref_2 + 1;
    elseif data( i , 5 ) == 4
        ref_4 = ref_4 + 1;
    elseif data( i , 5 ) == 6 
        ref_6 = ref_6 + 1;
    else
        disp('UNKNOWN REF NUMBER');
    end
end


% ZORADENIE A VYKRESLENIE  %
if ref_4 == 0 && ref_6 == 0
    
    [data_plus,data_minus,n_plus,n_minus] = pivot_split_counter(data);
    data_plus = bubble_sort( data_plus , n_plus );
    data_minus = bubble_sort( data_minus , n_minus );
    data_plus = flip(data_plus); % prevrati poradie riadkov
    
    data_sorted = [ data_minus ; data_plus ];
     
elseif ref_6 == 0
    [data_plus,data_minus,n_plus,n_minus] = manual_pivot_split_counter(data,line_points_1);
    
    [data_minus_r , data_minus_l , n_minus_r , n_minus_l] = pivot_split_counter(data_minus);
    data_minus_r = bubble_sort( data_minus_r , n_minus_r );
    data_minus_l = bubble_sort( data_minus_l , n_minus_l );
    data_minus_r = flip(data_minus_r); % prevrati poradie riadkov
    
    [data_plus_r , data_plus_l , n_plus_r , n_plus_l] = pivot_split_counter(data_plus);
    data_plus_r = bubble_sort( data_plus_r , n_plus_r );
    data_plus_l = bubble_sort( data_plus_l , n_plus_l );
    data_plus_r = flip(data_plus_r); % prevrati poradie riadkov
    
    data_sorted = [data_minus_l ; data_minus_r ; data_plus_l ; data_plus_r];
    
    line_1_active = 1;
     
else
    [data_plus,data_left,n_plus,n_left] = manual_pivot_split_counter(data,line_points_1);
    
    [data_l_plus,data_l_minus,n_l_plus,n_l_minus] = pivot_split_counter(data_left);
    data_l_plus = bubble_sort( data_l_plus , n_l_plus );
    data_l_minus = bubble_sort( data_l_minus , n_l_minus );
    data_l_plus = flip(data_l_plus); % prevrati poradie riadkov
    
    [data_right,data_mid,n_mid,n_right] = manual_pivot_split_counter(data_plus,line_points_2);
    [data_r_plus,data_r_minus,n_m_plus,n_m_minus] = pivot_split_counter(data_right);
    data_r_plus = bubble_sort( data_r_plus , n_m_plus );
    data_r_minus = bubble_sort( data_r_minus , n_m_minus );
    data_r_plus = flip(data_r_plus); % prevrati poradie riadkov
    
    [data_m_plus,data_m_minus,n_r_plus,n_r_minus] = pivot_split_counter(data_mid);
    data_m_plus = bubble_sort( data_m_plus , n_r_plus );
    data_m_minus = bubble_sort( data_m_minus , n_r_minus );
    data_m_plus = flip(data_m_plus); % prevrati poradie riadkov
    
    data_sorted = [data_l_minus ; data_l_plus ; data_m_minus ; data_m_plus; data_r_minus ; data_r_plus];
   
    line_1_active = 1;
    line_2_active = 1;
end

% VYKRESLENIE %

if line_1_active == plot_line && line_2_active == plot_line
    plot( line_points_1( : , 1 ),line_points_1( : , 2 ),'b*-',...
        line_points_2( : , 1 ),line_points_2( : , 2 ),'m*-');
    hold on;
elseif line_1_active == plot_line 
    plot( line_points_1( : , 1 ),line_points_1( : , 2 ),'bo-');
    hold on;   
end

ON_x = data_sorted( : , 1 );
ON_y = data_sorted( : , 2 );
OFF_x = data_sorted( : , 3 );
OFF_y = data_sorted( : , 4 );

plot(ON_x,ON_y,'g*-',OFF_x,OFF_y,'r*-');
xlabel('x [ mm ]');
ylabel('y [ mm ]');

if line_1_active == plot_line && line_2_active == plot_line
    legend('deliaca priamka 1','deliaca priamka 2','ON','OFF');
elseif line_1_active == plot_line 
    legend('deliaca priamka','ON','OFF');
else
    legend('ON','OFF');
end

   
function [d,d_ref_left] = position_detection( line_points , Px , Py )

global data;

Ax = line_points( 1 , 1 );
Ay = line_points( 1 , 2 );
Bx = line_points( 2 , 1 );
By = line_points( 2 , 2 );

d = (Px - Ax) * (By - Ay) - (Py - Ay) * (Bx - Ax); % hodnota pre urcenie polohy

min_x = min( data( : , 1 ) ); % bod ktory je najviac v lavo
d_ref_left = (min_x - Ax) * (By - Ay) - (0 - Ay) * (Bx - Ax); % referencna 
                                                   %hodnota pre lavu stranu

end

function [data_plus,data_minus,n_plus,n_minus] = manual_pivot_split_counter(data_in,line_points)

n_plus = 0; % pocitadlo pre body na lavo od priamky
n_minus = 0; % pocitadlo pre body na pravo od priamky
length_ = length( data_in( : , 1 ) ); % dlzka pola

[~,d_ref_left] = position_detection( line_points , 0 , 0 ); % referencna hodnota pre lave body

if d_ref_left < 0 % znamienka pre triedenie su nastavene pre zapornu hodnotu d_ref_left
    compensation = 1;
elseif d_red_left > 0 % pri kladnej hodnote d_ref_left sa kompenzuje prensasobenim -1
    compensation = -1;
else
    compensation = 1;
end
    

% spocita prvky na lavej a pravej strane od deliacej priamky %
for i = 1 : 1 : length_
    
    [d,~] = position_detection( line_points , data_in(i,1) , data_in(i,2) ); 
    % hodnota pre urcenie polohy
    
    if compensation * d > 0 % prvok je na pravej strane
        n_plus = n_plus + 1;
    elseif compensation * d < 0 % prvok je na lavej strane
        n_minus = n_minus + 1;
    else  % d == 0
        disp('PIVOT ERROR');
    end
end

% rozdeli prvky podla deliacej priamky %
data_plus = zeros( n_plus , 4); % alokacia pamate
data_minus = zeros( n_minus , 4); % alokacia pamate

for h = 1 : 2 : 3 % ON suradnice su v 1. a 2. stlpci, OFF v 3. a 4. stlpci
    
    data_plus_count = 1;
    data_minus_count = 1;

    for i = 1 : 1 : length_
        
        [d,~] = position_detection( line_points , data_in( i , h ) , data_in( i , h + 1 ) );
        
        if compensation * d > 0 % prvok je na pravej strane
            data_plus(data_plus_count , h : h + 1 ) = data_in( i , h : h + 1 );
            data_plus_count = data_plus_count + 1;
        else % compensation * d < 0 % prvok je na lavej strane
            data_minus(data_minus_count , h : h + 1 ) = data_in( i , h : h + 1 );
            data_minus_count = data_minus_count + 1;
        end
    end

end
end

function [data_plus,data_minus,n_plus,n_minus] = pivot_split_counter(data_in)
% pivot %
max_data = max( data_in( : , 2 ) ); % max hodnota y 
max_= zeros(2); % alokacia pamate
max_count = 1; % pocitadlo
length_ = length( data_in( : , 1 ) ); % dlzka vstupnych dat

for i = 1 : 1 : length_ % najde x polohy pre maxima y
    if max_data == data_in(i,2)
        max_(max_count) = data_in(i,1);
        max_count = max_count + 1;
    end
end
pivot = ( max_(1) + max_(2) ) / 2; % definuje pivot na x osi
n_plus = 0; % pocitadlo pre body na lavo od priamky
n_minus = 0; % pocitadlo pre body na pravo od priamky

% spocita prvky na lavej a pravej strane od pivot %
for i = 1 : 1 : length_
    if data_in(i,1) > pivot  % prvok je na pravej strane
        n_plus = n_plus + 1;
    elseif data_in(i,1) < pivot % prvok je na lavej strane
        n_minus = n_minus + 1;
    else
        disp('PIVOT ERROR');
    end
end

% rozdeli prvky podla pivot %
data_plus = zeros( n_plus , 4);
data_minus = zeros( n_minus , 4);

for h = 1 : 2 : 3 
    
    data_plus_count = 1;
    data_minus_count = 1;

    for i = 1 : 1 : length_
        if data_in( i , h ) > pivot % prvok je na pravej strane
            data_plus(data_plus_count , h : h + 1 ) = data_in( i , h : h + 1 );
            data_plus_count = data_plus_count + 1;
        else % data_in( i , h ) < pivot % prvok je na lavej strane
            data_minus(data_minus_count , h : h + 1 ) = data_in( i , h : h + 1 );
            data_minus_count = data_minus_count + 1;
        end
    end

end

end

function [data_in] = bubble_sort(data_in,length_)

for h = 2 : 2 : 4 % ON suradnice su v 1. a 2. stlpci, OFF v 3. a 4. stlpci
    
stop = 1;

while stop > 0
    stop = 0;
    for i = 1 : 1 : length_ - 1
        if data_in( i , h ) > data_in( i + 1 , h)
            temp = data_in( i + 1 , h - 1 : h );
            data_in( i + 1 , h - 1 : h ) = data_in( i , h - 1 : h );
            data_in( i ,  h - 1 : h ) = temp;
            stop = stop + 1;
        end
    end
end

end
[~,data_in] = zero_line_del(data_in,length_);
end

function  [length_,data_in] = zero_line_del(data_in,length_)

i = 1; % pocitadlo

while i <= length_ % mazanie riadkov obsahujucich iba nuly
    if(data_in(i,1)~=0 && data_in (i,2)~=0 && data_in(i,3)~=0 && data_in(i,4)~=0)
        i = i + 1;  
    else
        data_in(i,:)=[]; %odstranenie nuloveho riadku
    end 
    [length_ , ~] = size(data_in); %zmeranie poctu riadkov
end

end

function add_reference_number()
global data;
length_ = length(data(:,1));
data(:,5) = zeros(length_,1); % alokovanie pamate
%start = data(1,2); % hodnota start pre spustenie while loop

% for cyklus pridava referenciu do 5 stlpca data, na zaklade referencie sa 
% bude urcovat tvar pozadovanej krivky

for i = 1 : 1 : length_
    
    temp = data( i , 2 );
    same_y_count = 0;
    
    for j = 1 : 1 : length_
        
        if data( j , 2 ) == temp 
            same_y_count = same_y_count + 1; 
        end
        
    end
    data( i , 5 ) = same_y_count;
end

end
