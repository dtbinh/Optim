clear all
close all

%%
x = 1:100;
y = 1:100;
f = zeros(100,100);

for i = 1:100
    for j = 1:100
        f(i,j) = min(x(i),y(j));
    end
end    

figure
hold on
for i in x
    for j in y
        plot3(x(i),y(j),flip(f(i,j))
    end
end
