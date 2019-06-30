regulus_spica_px = 477;%477.3
regulus_arcturus_px = 538;%538.1
spica_arcturus_px = 294;

regulus_spica = 0.943547;
regulus_arcturus = 1.04206;
spica_arcturus = 0.572276;

factor_tmp(1) = regulus_spica / regulus_spica_px;
factor_tmp(2) = regulus_arcturus / regulus_arcturus_px;
factor_tmp(3) = spica_arcturus / spica_arcturus_px;

factor = mean(factor_tmp);
dfactor = sqrt(sum((factor_tmp - factor).^2) / (numel(factor_tmp)*(numel(factor_tmp)-1)));


%--------------------------------------------------------------------------
img = imread('lunar-halo.jpg');

xc = 623;
yc = 380;
dc = 10;

luna_tmp = img(yc-dc:yc+dc, xc-dc:xc+dc, :);
counter = 0;
for i = xc-dc:xc+dc
    for j = yc-dc:yc+dc
        counter = counter +1;
        luna_A(counter) = sum(img(j,i,:));
        luna_x(counter) = i;
        luna_y(counter) = j;
        luna_mesh(j,i) = sum(img(j,i,:));
    end
end

m = mean(luna_A);

luna_tmp2 = luna_tmp;
baricenter_X = 0;
baricenter_Y = 0;
peso = 0;
counter = 0;
for i = xc-dc:xc+dc
    for j = yc-dc:yc+dc
        counter = counter +1;
        if luna_A(counter) > m
            luna_tmp2(i - (xc-dc) + 1,j - (yc-dc) + 1,:) = luna_tmp(i - (xc-dc) + 1,j - (yc-dc) + 1,:);
            baricenter_X = baricenter_X + i * (luna_A(counter) - m);
            baricenter_Y = baricenter_Y + j * (luna_A(counter) - m);
            peso = peso + luna_A(counter) - m;
        else
            luna_tmp2(i - (xc-dc) + 1,j - (yc-dc) + 1,:) = 0;
        end
    end
end

baricenter_X = baricenter_X / peso;
baricenter_Y = baricenter_Y / peso;
imshow(luna_tmp2);
%figure, plot3(luna_x, luna_y, luna_A,'d');
figure, surf(luna_mesh(yc-dc:yc+dc, xc-dc:xc+dc))

xc = baricenter_X;
yc = baricenter_Y;

%--------------------------------------------------------------------------

load('luna.mat')

for i = 1:numel(x_r)
    r_r(i) = sqrt((x_r(i) - xc)^2 + (y_r(i) - yc)^2);
    t_r(i) = atan2(y_r(i) - yc, x_r(i) - xc);
end
for i = 1:numel(x_g)
    r_g(i) = sqrt((x_g(i) - xc)^2 + (y_g(i) - yc)^2);
    t_g(i) = atan2(y_g(i) - yc, x_g(i) - xc);
end
for i = 1:numel(x_b)
    r_b(i) = sqrt((x_b(i) - xc)^2 + (y_b(i) - yc)^2);
    t_b(i) = atan2(y_b(i) - yc, x_b(i) - xc);
end

r_r_m = mean(r_r);
r_g_m = mean(r_g);
r_b_m = mean(r_b);

dr_r_m = sqrt(sum((r_r - r_r_m).^2) / (numel(r_r)*(numel(r_r)-1)));
dr_g_m = sqrt(sum((r_g - r_g_m).^2) / (numel(r_g)*(numel(r_g)-1)));
dr_b_m = sqrt(sum((r_b - r_b_m).^2) / (numel(r_b)*(numel(r_b)-1)));

figure,
polarplot(-t_r, r_r - r_r_m, '+r'); hold on;
polarplot(-t_g, r_g - r_g_m, '+g');
polarplot(-t_b, r_b - r_b_m, '+b');
polarplot(0:0.001:2*pi, (0:0.001:2*pi)*0, 'k');
rlim([-7 7])

n = @(r, f) sin((r + f)/2)/sin(f/2);
n_lati = 12;
phi = -pi + 2 * pi*(n_lati-2)/n_lati;
%phi = pi/2;
n_r = n(r_r_m * factor, phi);
n_g = n(r_g_m * factor, phi);
n_b = n(r_b_m * factor, phi);

dn_r = n((r_r_m + dr_r_m) * (factor + dfactor), 60/180*pi) - n_r;
dn_g = n((r_g_m + dr_r_m) * (factor + dfactor), 60/180*pi) - n_g;
dn_b = n((r_b_m + dr_r_m) * (factor + dfactor), 60/180*pi) - n_b;











