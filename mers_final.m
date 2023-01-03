clear; clf; format long
population = 50220000;
MERS_I = load('mers_infective.m'); MERS_R = load('mers_recovered.m'); 

date = datenum('15-jul-2015') - datenum('20-may-2015');

%% day1~day2
date1 = datenum('21-may-2015') - datenum('20-may-2015');
sc = 100; Nt = sc*date1; t = linspace(0, date1, Nt+1); dt = t(2)-t(1);

s1 = zeros(1,Nt); i1 = zeros(1,Nt); r1 = zeros(1,Nt);
s1(1) = population-2; i1(1) = 2; r1(1) = 0;

b = 0.0000001250; bf1 = 0.065; k1 = 0.0; 
for j = 1:Nt
    s1(j+1) = s1(j) - b*bf1*s1(j)*i1(j)*dt;
    i1(j+1) = i1(j) + (b*bf1*s1(j)*i1(j) - (k1)*i1(j))*dt;
    r1(j+1) = r1(j) + k1*i1(j)*dt;
end
figure(1); clf;
hold on
plot(t, i1, 'b-', t, r1, 'g-');
plot(0:date1, MERS_I(1:2),'ro-', 0:date1, MERS_R(1:2),'mo-');
legend('Infected','Recovered',2 );


%% day2~day3
clear Nt t dt b; format long;
date2 = datenum('22-may-2015') - datenum('21-may-2015');
Nt = sc*date2; t = linspace(0, date2, Nt+1); dt = t(2)-t(1);

s2 = zeros(1,Nt); i2 = zeros(1,Nt); r2 = zeros(1,Nt);
s2(1) = population-i1(end); i2(1) = i1(end); r2(1) = r1(end);

b = 0.0000001250; bf2 = 0.0; k2 = 0.0; 
for j = 1:Nt
    s2(j+1) = s2(j) - b*bf2*s2(j)*i2(j)*dt;
    i2(j+1) = i2(j) + (b*bf2*s2(j)*i2(j) - (k2)*i2(j))*dt;
    r2(j+1) = r2(j) + k2*i2(j)*dt;
end
figure(2); clf;
hold on
plot(t, i2, 'b-', t, r2, 'g-');
plot(0:date2, MERS_I(2:3),'ro-', 0:date2, MERS_R(2:3),'mo-');
legend('Infected','Recovered',2);

%% day3~day4
clear Nt t dt b; format long;
date3 = datenum('23-may-2015') - datenum('22-may-2015');
Nt = sc*date3; t = linspace(0, date3, Nt+1); dt = t(2)-t(1);

s3 = zeros(1,Nt); i3 = zeros(1,Nt); r3 = zeros(1,Nt);
s3(1) = population-i2(end); i3(1) = i2(end); r3(1) = r2(end);

b = 0.0000001250; bf3 = 0.0; k3 = 0.0; 
for j = 1:Nt
    s3(j+1) = s3(j) - b*bf3*s3(j)*i3(j)*dt;
    i3(j+1) = i3(j) + (b*bf3*s3(j)*i3(j) - (k3)*i3(j))*dt;
    r3(j+1) = r3(j) + k3*i3(j)*dt;
end
figure(3); clf;
hold on
plot(t, i3, 'b-', t, r3, 'g-');
plot(0:date3, MERS_I(3:4),'ro-', 0:date3, MERS_R(3:4),'mo-');
legend('Infected','Recovered',2);

%% day4~day5
clear Nt t dt b; format long;
date4 = datenum('24-may-2015') - datenum('23-may-2015');
Nt = sc*date4; t = linspace(0, date4, Nt+1); dt = t(2)-t(1);

s4 = zeros(1,Nt); i4 = zeros(1,Nt); r4 = zeros(1,Nt);
s4(1) = population-i3(end); i4(1) = i3(end); r4(1) = r3(end);

b = 0.0000001250; bf4 = 0.0; k4 = 0.0; 
for j = 1:Nt
    s4(j+1) = s4(j) - b*bf4*s4(j)*i4(j)*dt;
    i4(j+1) = i4(j) + (b*bf4*s4(j)*i4(j) - (k4)*i4(j))*dt;
    r4(j+1) = r4(j) + k4*i4(j)*dt;
end
figure(4); clf;
hold on
plot(t, i4, 'b-', t, r4, 'g-');
plot(0:date4, MERS_I(4:5),'ro-', 0:date4, MERS_R(4:5),'mo-');
legend('Infected','Recovered',2);


%% day5~day6
clear Nt t dt b; format long;
date5 = datenum('25-may-2015') - datenum('24-may-2015');
Nt = sc*date5; t = linspace(0, date5, Nt+1); dt = t(2)-t(1);

s5 = zeros(1,Nt); i5 = zeros(1,Nt); r5 = zeros(1,Nt);
s5(1) = population-i4(end); i5(1) = i4(end); r5(1) = r4(end);

b = 0.0000001250; bf5 = 0.0; k5 = 0.0; 
for j = 1:Nt
    s5(j+1) = s5(j) - b*bf5*s5(j)*i5(j)*dt;
    i5(j+1) = i5(j) + (b*bf5*s5(j)*i5(j) - (k5)*i5(j))*dt;
    r5(j+1) = r5(j) + k5*i5(j)*dt;
end
figure(5); clf;
hold on
plot(t, i5, 'b-', t, r5, 'g-');
plot(0:date5, MERS_I(5:6),'ro-', 0:date5, MERS_R(5:6),'mo-');
legend('Infected','Recovered',2);

%% day6~day7
clear Nt t dt b; format long;
date6 = datenum('26-may-2015') - datenum('25-may-2015');
Nt = sc*date6; t = linspace(0, date6, Nt+1); dt = t(2)-t(1);

s6 = zeros(1,Nt); i6 = zeros(1,Nt); r6 = zeros(1,Nt);
s6(1) = population-i5(end); i6(1) = i5(end); r6(1) = r5(end);

b = 0.0000001250; bf6 = 0.081; k6 = 0.0; 
for j = 1:Nt
    s6(j+1) = s6(j) - b*bf6*s6(j)*i6(j)*dt;
    i6(j+1) = i6(j) + (b*bf6*s6(j)*i6(j) - (k6)*i6(j))*dt;
    r6(j+1) = r6(j) + k6*i6(j)*dt;
end
figure(6); clf;
hold on
plot(t, i6, 'b-', t, r6, 'g-');
plot(0:date6, MERS_I(6:7),'ro-', 0:date6, MERS_R(6:7),'mo-');
legend('Infected','Recovered',2);


%% day7~day8
clear Nt t dt b; format long;
date7 = datenum('27-may-2015') - datenum('26-may-2015');
Nt = sc*date7; t = linspace(0, date7, Nt+1); dt = t(2)-t(1);

s7 = zeros(1,Nt); i7 = zeros(1,Nt); r7 = zeros(1,Nt);
s7(1) = population-i6(end); i7(1) = i6(end); r7(1) = r6(end);

b = 0.0000001250; bf7 = 0.00; k7 = 0.0; 
for j = 1:Nt
    s7(j+1) = s7(j) - b*bf7*s7(j)*i7(j)*dt;
    i7(j+1) = i7(j) + (b*bf7*s7(j)*i7(j) - (k7)*i7(j))*dt;
    r7(j+1) = r7(j) + k7*i7(j)*dt;
end
figure(7); clf;
hold on
plot(t, i7, 'b-', t, r7, 'g-');
plot(0:date7, MERS_I(7:8),'ro-', 0:date7, MERS_R(7:8),'mo-');
legend('Infected','Recovered',2);

%% day8~day9
clear Nt t dt b; format long;
date8 = datenum('28-may-2015') - datenum('27-may-2015');
Nt = sc*date8; t = linspace(0, date8, Nt+1); dt = t(2)-t(1);

s8 = zeros(1,Nt); i8 = zeros(1,Nt); r8 = zeros(1,Nt);
s8(1) = population-i7(end); i8(1) = i7(end); r8(1) = r7(end);

b = 0.0000001250; bf8 = 0.054; k8 = 0.0; 
for j = 1:Nt
    s8(j+1) = s8(j) - b*bf8*s8(j)*i8(j)*dt;
    i8(j+1) = i8(j) + (b*bf8*s8(j)*i8(j) - (k8)*i8(j))*dt;
    r8(j+1) = r8(j) + k8*i8(j)*dt;
end
figure(8); clf;
hold on
plot(t, i8, 'b-', t, r8, 'g-');
plot(0:date8, MERS_I(8:9),'ro-', 0:date8, MERS_R(8:9),'mo-');
legend('Infected','Recovered',2);


%% day9~day10
clear Nt t dt b; format long;
date9 = datenum('29-may-2015') - datenum('28-may-2015');
Nt = sc*date9; t = linspace(0, date9, Nt+1); dt = t(2)-t(1);

s9 = zeros(1,Nt); i9 = zeros(1,Nt); r9 = zeros(1,Nt);
s9(1) = population-i8(end); i9(1) = i8(end); r9(1) = r8(end);

b = 0.0000001250; bf9 = 0.099; k9 = 0.0; 
for j = 1:Nt
    s9(j+1) = s9(j) - b*bf9*s9(j)*i9(j)*dt;
    i9(j+1) = i9(j) + (b*bf9*s9(j)*i9(j) - (k9)*i9(j))*dt;
    r9(j+1) = r9(j) + k9*i9(j)*dt;
end
figure(9); clf;
hold on
plot(t, i9, 'b-', t, r9, 'g-');
plot(0:date9, MERS_I(9:10),'ro-', 0:date9, MERS_R(9:10),'mo-');
legend('Infected','Recovered',2);

%% day10~day11
clear Nt t dt b; format long;
date10 = datenum('30-may-2015') - datenum('29-may-2015');
Nt = sc*date10; t = linspace(0, date10, Nt+1); dt = t(2)-t(1);

s10 = zeros(1,Nt); i10 = zeros(1,Nt); r10 = zeros(1,Nt);
s10(1) = population-i9(end); i10(1) = i9(end); r10(1) = r9(end);

b = 0.0000001250; bf10 = 0.023; k10 = 0.0; 
for j = 1:Nt
    s10(j+1) = s10(j) - b*bf10*s10(j)*i10(j)*dt;
    i10(j+1) = i10(j) + (b*bf10*s10(j)*i10(j) - (k10)*i10(j))*dt;
    r10(j+1) = r10(j) + k10*i10(j)*dt;
end
figure(10); clf;
hold on
plot(t, i10, 'b-', t, r10, 'g-');
plot(0:date10, MERS_I(10:11),'ro-', 0:date10, MERS_R(10:11),'mo-');
legend('Infected','Recovered',2);


%% day11~day12
clear Nt t dt b; format long;
date11 = datenum('31-may-2015') - datenum('30-may-2015');
Nt = sc*date11; t = linspace(0, date11, Nt+1); dt = t(2)-t(1);

s11 = zeros(1,Nt); i11 = zeros(1,Nt); r11 = zeros(1,Nt);
s11(1) = population-i10(end); i11(1) = i10(end); r11(1) = r10(end);

b = 0.0000001250; bf11 = 0.028; k11 = 0.0; 
for j = 1:Nt
    s11(j+1) = s11(j) - b*bf11*s11(j)*i11(j)*dt;
    i11(j+1) = i11(j) + (b*bf11*s11(j)*i11(j) - (k11)*i11(j))*dt;
    r11(j+1) = r11(j) + k11*i11(j)*dt;
end
figure(11); clf;
hold on
plot(t, i11, 'b-', t, r11, 'g-');
plot(0:date11, MERS_I(11:12),'ro-', 0:date11, MERS_R(11:12),'mo-');
legend('Infected','Recovered',2);

%% day12~day13
clear Nt t dt b; format long;
date12 = datenum('1-june-2015') - datenum('31-may-2015');
Nt = sc*date12; t = linspace(0, date12, Nt+1); dt = t(2)-t(1);

s12 = zeros(1,Nt); i12 = zeros(1,Nt); r12 = zeros(1,Nt);
s12(1) = population-i11(end); i12(1) = i11(end); r12(1) = r11(end);

b = 0.0000001250; bf12 = 0.055; k12 = 0.05; 
for j = 1:Nt
    s12(j+1) = s12(j) - b*bf12*s12(j)*i12(j)*dt;
    i12(j+1) = i12(j) + (b*bf12*s12(j)*i12(j) - (k12)*i12(j))*dt;
    r12(j+1) = r12(j) + k12*i12(j)*dt;
end
figure(12); clf;
hold on
plot(t, i12, 'b-', t, r12, 'g-');
plot(0:date12, MERS_I(12:13),'ro-', 0:date12, MERS_R(12:13),'mo-');
legend('Infected','Recovered',2);

%% day13~day14
clear Nt t dt b; format long;
date13 = datenum('2-june-2015') - datenum('1-june-2015');
Nt = sc*date13; t = linspace(0, date13, Nt+1); dt = t(2)-t(1);

s13 = zeros(1,Nt); i13 = zeros(1,Nt); r13 = zeros(1,Nt);
s13(1) = population-i12(end); i13(1) = i12(end); r13(1) = r12(end);

b = 0.0000001250; bf13 = 0.03; k13 = 0.0; 
for j = 1:Nt
    s13(j+1) = s13(j) - b*bf13*s13(j)*i13(j)*dt;
    i13(j+1) = i13(j) + (b*bf13*s13(j)*i13(j) - (k13)*i13(j))*dt;
    r13(j+1) = r13(j) + k13*i13(j)*dt;
end
figure(13); clf;
hold on
plot(t, i13, 'b-', t, r13, 'g-');
plot(0:date13, MERS_I(13:14),'ro-', 0:date13, MERS_R(13:14),'mo-');
legend('Infected','Recovered',2);


%% day14~day15
clear Nt t dt b; format long;
date14 = datenum('3-june-2015') - datenum('2-june-2015');
Nt = sc*date14; t = linspace(0, date14, Nt+1); dt = t(2)-t(1);

s14 = zeros(1,Nt); i14 = zeros(1,Nt); r14 = zeros(1,Nt);
s14(1) = population-i13(end); i14(1) = i13(end); r14(1) = r13(end);

b = 0.0000001250; bf14 = 0.0008; k14 = 0.07; 
for j = 1:Nt
    s14(j+1) = s14(j) - b*bf14*s13(j)*i14(j)*dt;
    i14(j+1) = i14(j) + (b*bf14*s13(j)*i14(j) - (k14)*i14(j))*dt;
    r14(j+1) = r14(j) + k14*i14(j)*dt;
end
figure(14); clf;
hold on
plot(t, i14, 'b-', t, r14, 'g-');
plot(0:date14, MERS_I(14:15),'ro-', 0:date14, MERS_R(14:15),'mo-');
legend('Infected','Recovered',2);


%% day15~day16
clear Nt t dt b; format long;
date15 = datenum('4-june-2015') - datenum('3-june-2015');
Nt = sc*date15; t = linspace(0, date15, Nt+1); dt = t(2)-t(1);

s15 = zeros(1,Nt); i15 = zeros(1,Nt); r15 = zeros(1,Nt);
s15(1) = population-i14(end); i15(1) = i14(end); r15(1) = r14(end);

b = 0.0000001250; bf15 = 0.032; k15 = 0.04; 
for j = 1:Nt
    s15(j+1) = s15(j) - b*bf15*s15(j)*i15(j)*dt;
    i15(j+1) = i15(j) + (b*bf15*s15(j)*i15(j) - (k15)*i15(j))*dt;
    r15(j+1) = r15(j) + k15*i15(j)*dt;
end
figure(15); clf;
hold on
plot(t, i15, 'b-', t, r15, 'g-');
plot(0:date15, MERS_I(15:16),'ro-', 0:date15, MERS_R(15:16),'mo-');
legend('Infected','Recovered',2);


%% day16~day17
clear Nt t dt b; format long;
date16 = datenum('5-june-2015') - datenum('4-june-2015');
Nt = sc*date16; t = linspace(0, date16, Nt+1); dt = t(2)-t(1);

s16 = zeros(1,Nt); i16 = zeros(1,Nt); r16 = zeros(1,Nt);
s16(1) = population-i15(end); i16(1) = i15(end); r16(1) = r15(end);

b = 0.0000001250; bf16 = 0.03; k16 = 0.06; 
for j = 1:Nt
    s16(j+1) = s16(j) - b*bf16*s16(j)*i16(j)*dt;
    i16(j+1) = i16(j) + (b*bf16*s16(j)*i16(j) - (k16)*i16(j))*dt;
    r16(j+1) = r16(j) + k16*i16(j)*dt;
end
figure(16); clf;
hold on
plot(t, i16, 'b-', t, r16, 'g-');
plot(0:date16, MERS_I(16:17),'ro-', 0:date16, MERS_R(16:17),'mo-');
legend('Infected','Recovered',2);


%% day17~day18
clear Nt t dt b; format long;
date17 = datenum('6-june-2015') - datenum('5-june-2015');
Nt = sc*date17; t = linspace(0, date17, Nt+1); dt = t(2)-t(1);

s17 = zeros(1,Nt); i17 = zeros(1,Nt); r17 = zeros(1,Nt);
s17(1) = population-i16(end); i17(1) = i16(end); r17(1) = r16(end);

b = 0.0000001250; bf17 = 0.075; k17 = 0.0; 
for j = 1:Nt
    s17(j+1) = s17(j) - b*bf17*s17(j)*i17(j)*dt;
    i17(j+1) = i17(j) + (b*bf17*s17(j)*i17(j) - (k17)*i17(j))*dt;
    r17(j+1) = r17(j) + k17*i17(j)*dt;
end
figure(17); clf;
hold on
plot(t, i17, 'b-', t, r17, 'g-');
plot(0:date17, MERS_I(17:18),'ro-', 0:date17, MERS_R(17:18),'mo-');
legend('Infected','Recovered',2);

%% day18~day19
clear Nt t dt b; format long;
date18 = datenum('7-june-2015') - datenum('6-june-2015');
Nt = sc*date18; t = linspace(0, date18, Nt+1); dt = t(2)-t(1);

s18 = zeros(1,Nt); i18 = zeros(1,Nt); r18 = zeros(1,Nt);
s18(1) = population-i17(end); i18(1) = i17(end); r18(1) = r17(end);

b = 0.0000001250; bf18 = 0.053; k18 = 0.0; 
for j = 1:Nt
    s18(j+1) = s18(j) - b*bf18*s18(j)*i18(j)*dt;
    i18(j+1) = i18(j) + (b*bf18*s18(j)*i18(j) - (k18)*i18(j))*dt;
    r18(j+1) = r18(j) + k18*i18(j)*dt;
end
figure(18); clf;
hold on
plot(t, i18, 'b-', t, r18, 'g-');
plot(0:date18, MERS_I(18:19),'ro-', 0:date18, MERS_R(18:19),'mo-');
legend('Infected','Recovered',2);


%% day19~day20
clear Nt t dt b; format long;
date19 = datenum('8-june-2015') - datenum('7-june-2015');
Nt = sc*date19; t = linspace(0, date19, Nt+1); dt = t(2)-t(1);

s19 = zeros(1,Nt); i19 = zeros(1,Nt); r19 = zeros(1,Nt);
s19(1) = population-i18(end); i19(1) = i18(end); r19(1) = r18(end);

b = 0.0000001250; bf19 = 0.015; k19 = 0.03; 
for j = 1:Nt
    s19(j+1) = s19(j) - b*bf19*s19(j)*i19(j)*dt;
    i19(j+1) = i19(j) + (b*bf19*s19(j)*i19(j) - (k19)*i19(j))*dt;
    r19(j+1) = r19(j) + k19*i19(j)*dt;
end
figure(19); clf;
hold on
plot(t, i19, 'b-', t, r19, 'g-');
plot(0:date19, MERS_I(19:20),'ro-', 0:date19, MERS_R(19:20),'mo-');
legend('Infected','Recovered',2);


%% day20~day21
clear Nt t dt b; format long;
date20 = datenum('9-june-2015') - datenum('8-june-2015');
Nt = sc*date20; t = linspace(0, date20, Nt+1); dt = t(2)-t(1);

s20 = zeros(1,Nt); i20 = zeros(1,Nt); r20 = zeros(1,Nt);
s20(1) = population-i19(end); i20(1) = i19(end); r20(1) = r19(end);

b = 0.0000001250; bf20 = 0.023; k20 = 0.015; 
for j = 1:Nt
    s20(j+1) = s20(j) - b*bf20*s20(j)*i20(j)*dt;
    i20(j+1) = i20(j) + (b*bf20*s20(j)*i20(j) - (k20)*i20(j))*dt;
    r20(j+1) = r20(j) + k20*i20(j)*dt;
end
figure(20); clf;
hold on
plot(t, i20, 'b-', t, r20, 'g-');
plot(0:date20, MERS_I(20:21),'ro-', 0:date20, MERS_R(20:21),'mo-');
legend('Infected','Recovered',2);


%% day21~day22
clear Nt t dt b; format long;
date21 = datenum('10-june-2015') - datenum('9-june-2015');
Nt = sc*date21; t = linspace(0, date21, Nt+1); dt = t(2)-t(1);

s21 = zeros(1,Nt); i21 = zeros(1,Nt); r21 = zeros(1,Nt);
s21(1) = population-i20(end); i21(1) = i20(end); r21(1) = r20(end);

b = 0.0000001250; bf21 = 0.02; k21 = 0.025; 
for j = 1:Nt
    s21(j+1) = s21(j) - b*bf21*s21(j)*i21(j)*dt;
    i21(j+1) = i21(j) + (b*bf21*s21(j)*i21(j) - (k21)*i21(j))*dt;
    r21(j+1) = r21(j) + k21*i21(j)*dt;
end
figure(21); clf;
hold on
plot(t, i21, 'b-', t, r21, 'g-');
plot(0:date21, MERS_I(21:22),'ro-', 0:date21, MERS_R(21:22),'mo-');
legend('Infected','Recovered',2);

%% day22~day23
clear Nt t dt b; format long;
date22 = datenum('11-june-2015') - datenum('10-june-2015');
Nt = sc*date22; t = linspace(0, date22, Nt+1); dt = t(2)-t(1);

s22 = zeros(1,Nt); i22 = zeros(1,Nt); r22 = zeros(1,Nt);
s22(1) = population-i21(end); i22(1) = i21(end); r22(1) = r21(end);

b = 0.0000001250; bf22 = 0.005; k22 = 0.035; 
for j = 1:Nt
    s22(j+1) = s22(j) - b*bf22*s22(j)*i22(j)*dt;
    i22(j+1) = i22(j) + (b*bf22*s22(j)*i22(j) - (k22)*i22(j))*dt;
    r22(j+1) = r22(j) + k22*i22(j)*dt;
end
figure(22); clf;
hold on
plot(t, i22, 'b-', t, r22, 'g-');
plot(0:date22, MERS_I(22:23),'ro-', 0:date22, MERS_R(22:23),'mo-');
legend('Infected','Recovered',2);


%% day23~day24
clear Nt t dt b; format long;
date23 = datenum('12-june-2015') - datenum('11-june-2015');
Nt = sc*date23; t = linspace(0, date23, Nt+1); dt = t(2)-t(1);

s23 = zeros(1,Nt); i23 = zeros(1,Nt); r23 = zeros(1,Nt);
s23(1) = population-i22(end); i23(1) = i22(end); r23(1) = r22(end);

b = 0.0000001250; bf23 = 0.017; k23 = 0.043; 
for j = 1:Nt
    s23(j+1) = s23(j) - b*bf23*s23(j)*i23(j)*dt;
    i23(j+1) = i23(j) + (b*bf23*s23(j)*i23(j) - (k23)*i23(j))*dt;
    r23(j+1) = r23(j) + k23*i23(j)*dt;
end
figure(23); clf;
hold on
plot(t, i23, 'b-', t, r23, 'g-');
plot(0:date23, MERS_I(23:24),'ro-', 0:date23, MERS_R(23:24),'mo-');
legend('Infected','Recovered',2);


%% day24~day25
clear Nt t dt b; format long;
date24 = datenum('13-june-2015') - datenum('12-june-2015');
Nt = sc*date24; t = linspace(0, date24, Nt+1); dt = t(2)-t(1);

s24 = zeros(1,Nt); i24 = zeros(1,Nt); r24 = zeros(1,Nt);
s24(1) = population-i23(end); i24(1) = i23(end); r24(1) = r23(end);

b = 0.0000001250; bf24 = 0.01; k24 = 0.023; 
for j = 1:Nt
    s24(j+1) = s24(j) - b*bf24*s24(j)*i24(j)*dt;
    i24(j+1) = i24(j) + (b*bf24*s24(j)*i24(j) - (k24)*i24(j))*dt;
    r24(j+1) = r24(j) + k24*i24(j)*dt;
end
figure(24); clf;
hold on
plot(t, i24, 'b-', t, r24, 'g-');
plot(0:date24, MERS_I(24:25),'ro-', 0:date24, MERS_R(24:25),'mo-');
legend('Infected','Recovered',2);

%% day25~day26
clear Nt t dt b; format long;
date25 = datenum('14-june-2015') - datenum('13-june-2015');
Nt = sc*date25; t = linspace(0, date25, Nt+1); dt = t(2)-t(1);

s25 = zeros(1,Nt); i25 = zeros(1,Nt); r25 = zeros(1,Nt);
s25(1) = population-i24(end); i25(1) = i24(end); r25(1) = r24(end);

b = 0.0000001250; bf25 = 0.008; k25 = 0.05; 
for j = 1:Nt
    s25(j+1) = s25(j) - b*bf25*s25(j)*i25(j)*dt;
    i25(j+1) = i25(j) + (b*bf25*s25(j)*i25(j) - (k25)*i25(j))*dt;
    r25(j+1) = r25(j) + k25*i25(j)*dt;
end
figure(25); clf;
hold on
plot(t, i25, 'b-', t, r25, 'g-');
plot(0:date25, MERS_I(25:26),'ro-', 0:date25, MERS_R(25:26),'mo-');
legend('Infected','Recovered',2);


%% day26~day27
clear Nt t dt b; format long;
date26 = datenum('15-june-2015') - datenum('14-june-2015');
Nt = sc*date26; t = linspace(0, date26, Nt+1); dt = t(2)-t(1);

s26 = zeros(1,Nt); i26 = zeros(1,Nt); r26 = zeros(1,Nt);
s26(1) = population-i25(end); i26(1) = i25(end); r26(1) = r25(end);

b = 0.0000001250; bf26 = 0.006; k26 = 0.053; 
for j = 1:Nt
    s26(j+1) = s26(j) - b*bf26*s26(j)*i26(j)*dt;
    i26(j+1) = i26(j) + (b*bf26*s26(j)*i26(j) - (k26)*i26(j))*dt;
    r26(j+1) = r26(j) + k26*i26(j)*dt;
end
figure(26); clf;
hold on
plot(t, i26, 'b-', t, r26, 'g-');
plot(0:date26, MERS_I(26:27),'ro-', 0:date26, MERS_R(26:27),'mo-');
legend('Infected','Recovered',2);

%% day27~day28
clear Nt t dt b; format long;
date27 = datenum('16-june-2015') - datenum('15-june-2015');
Nt = sc*date27; t = linspace(0, date27, Nt+1); dt = t(2)-t(1);

s27 = zeros(1,Nt); i27 = zeros(1,Nt); r27 = zeros(1,Nt);
s27(1) = population-i26(end); i27(1) = i26(end); r27(1) = r26(end);

b = 0.0000001250; bf27 = 0.01; k27 = 0.015; 
for j = 1:Nt
    s27(j+1) = s27(j) - b*bf27*s27(j)*i27(j)*dt;
    i27(j+1) = i27(j) + (b*bf27*s27(j)*i27(j) - (k27)*i27(j))*dt;
    r27(j+1) = r27(j) + k27*i27(j)*dt;
end
figure(27); clf;
hold on
plot(t, i27, 'b-', t, r27, 'g-');
plot(0:date27, MERS_I(27:28),'ro-', 0:date27, MERS_R(27:28),'mo-');
legend('Infected','Recovered',2);

%% day28~day29
clear Nt t dt b; format long;
date28 = datenum('17-june-2015') - datenum('16-june-2015');
Nt = sc*date28; t = linspace(0, date28, Nt+1); dt = t(2)-t(1);

s28 = zeros(1,Nt); i28 = zeros(1,Nt); r28 = zeros(1,Nt);
s28(1) = population-i27(end); i28(1) = i27(end); r28(1) = r27(end);

b = 0.0000001250; bf28 = 0.0042; k28 = 0.075; 
for j = 1:Nt
    s28(j+1) = s28(j) - b*bf28*s28(j)*i28(j)*dt;
    i28(j+1) = i28(j) + (b*bf28*s28(j)*i28(j) - (k28)*i28(j))*dt;
    r28(j+1) = r28(j) + k28*i28(j)*dt;
end
figure(28); clf;
hold on
plot(t, i28, 'b-', t, r28, 'g-');
plot(0:date28, MERS_I(28:29),'ro-', 0:date28, MERS_R(28:29),'mo-');
legend('Infected','Recovered',2);

%% day29~day30
clear Nt t dt b; format long;
date29 = datenum('18-june-2015') - datenum('17-june-2015');
Nt = sc*date29; t = linspace(0, date29, Nt+1); dt = t(2)-t(1);

s29 = zeros(1,Nt); i29 = zeros(1,Nt); r29 = zeros(1,Nt);
s29(1) = population-i28(end); i29(1) = i28(end); r29(1) = r28(end);

b = 0.0000001250; bf29 = 0.0004; k29 = 0.06; 
for j = 1:Nt
    s29(j+1) = s29(j) - b*bf29*s29(j)*i29(j)*dt;
    i29(j+1) = i29(j) + (b*bf29*s29(j)*i29(j) - (k29)*i29(j))*dt;
    r29(j+1) = r29(j) + k29*i29(j)*dt;
end
figure(29); clf;
hold on
plot(t, i29, 'b-', t, r29, 'g-');
plot(0:date29, MERS_I(29:30),'ro-', 0:date29, MERS_R(29:30),'mo-');
legend('Infected','Recovered',2);

%% day30~day31
clear Nt t dt b; format long;
date30 = datenum('19-june-2015') - datenum('18-june-2015');
Nt = sc*date30; t = linspace(0, date30, Nt+1); dt = t(2)-t(1);

s30 = zeros(1,Nt); i30 = zeros(1,Nt); r30 = zeros(1,Nt);
s30(1) = population-i29(end); i30(1) = i29(end); r30(1) = r29(end);

b = 0.0000001250; bf30 = 0.0000001; k30 = 0.05; 
for j = 1:Nt
    s30(j+1) = s30(j) - b*bf30*s30(j)*i30(j)*dt;
    i30(j+1) = i30(j) + (b*bf30*s30(j)*i30(j) - (k30)*i30(j))*dt;
    r30(j+1) = r30(j) + k30*i30(j)*dt;
end
figure(30); clf;
hold on
plot(t, i30, 'b-', t, r30, 'g-');
plot(0:date30, MERS_I(30:31),'ro-', 0:date30, MERS_R(30:31),'mo-');
legend('Infected','Recovered',2);

%% day31~day32
clear Nt t dt b; format long;
date31 = datenum('20-june-2015') - datenum('19-june-2015');
Nt = sc*date31; t = linspace(0, date31, Nt+1); dt = t(2)-t(1);

s31 = zeros(1,Nt); i31 = zeros(1,Nt); r31 = zeros(1,Nt);
s31(1) = population-i22(end); i31(1) = i30(end); r31(1) = r30(end);

b = 0.0000001250; bf31 = 0.005; k31 = 0.08; 
for j = 1:Nt
    s31(j+1) = s31(j) - b*bf31*s31(j)*i31(j)*dt;
    i31(j+1) = i31(j) + (b*bf31*s31(j)*i31(j) - (k31)*i31(j))*dt;
    r31(j+1) = r31(j) + k31*i31(j)*dt;
end
figure(31); clf;
hold on
plot(t, i31, 'b-', t, r31, 'g-');
plot(0:date31, MERS_I(31:32),'ro-', 0:date31, MERS_R(31:32),'mo-');
legend('Infected','Recovered',2);

%% day32~day33
clear Nt t dt b; format long;
date32 = datenum('21-june-2015') - datenum('20-june-2015');
Nt = sc*date32; t = linspace(0, date32, Nt+1); dt = t(2)-t(1);

s32 = zeros(1,Nt); i32 = zeros(1,Nt); r32 = zeros(1,Nt);
s32(1) = population-i31(end); i32(1) = i31(end); r32(1) = r31(end);

b = 0.0000001250; bf32 = 0.005; k32 = 0.095; 
for j = 1:Nt
    s32(j+1) = s32(j) - b*bf32*s32(j)*i32(j)*dt;
    i32(j+1) = i32(j) + (b*bf32*s32(j)*i32(j) - (k32)*i32(j))*dt;
    r32(j+1) = r32(j) + k32*i32(j)*dt;
end
figure(32); clf;
hold on
plot(t, i32, 'b-', t, r32, 'g-');
plot(0:date32, MERS_I(32:33),'ro-', 0:date32, MERS_R(32:33),'mo-');
legend('Infected','Recovered',2);


%% day33~day34
clear Nt t dt b; format long;
date33 = datenum('22-june-2015') - datenum('21-june-2015');
Nt = sc*date33; t = linspace(0, date33, Nt+1); dt = t(2)-t(1);

s33 = zeros(1,Nt); i33 = zeros(1,Nt); r33 = zeros(1,Nt);
s33(1) = population-i32(end); i33(1) = i32(end); r33(1) = r32(end);

b = 0.0000001250; bf33 = 0.005; k33 = 0.04; 
for j = 1:Nt
    s33(j+1) = s33(j) - b*bf33*s33(j)*i33(j)*dt;
    i33(j+1) = i33(j) + (b*bf33*s33(j)*i33(j) - (k33)*i33(j))*dt;
    r33(j+1) = r33(j) + k33*i33(j)*dt;
end
figure(33); clf;
hold on
plot(t, i33, 'b-', t, r33, 'g-');
plot(0:date33, MERS_I(33:34),'ro-', 0:date33, MERS_R(33:34),'mo-');
legend('Infected','Recovered',2);


%% day34~day35
clear Nt t dt b; format long;
date34 = datenum('23-june-2015') - datenum('22-june-2015');
Nt = sc*date34; t = linspace(0, date34, Nt+1); dt = t(2)-t(1);

s34 = zeros(1,Nt); i34 = zeros(1,Nt); r34 = zeros(1,Nt);
s34(1) = population-i33(end); i34(1) = i33(end); r34(1) = r33(end);

b = 0.0000001250; bf34 = 0.0072; k34 = 0.145; 
for j = 1:Nt
    s34(j+1) = s34(j) - b*bf34*s34(j)*i34(j)*dt;
    i34(j+1) = i34(j) + (b*bf34*s34(j)*i34(j) - (k34)*i34(j))*dt;
    r34(j+1) = r34(j) + k34*i34(j)*dt;
end
figure(34); clf;
hold on
plot(t, i34, 'b-', t, r34, 'g-');
plot(0:date34, MERS_I(34:35),'ro-', 0:date34, MERS_R(34:35),'mo-');
legend('Infected','Recovered',2);

%% day35~day36
clear Nt t dt b; format long;
date35 = datenum('24-june-2015') - datenum('23-june-2015');
Nt = sc*date35; t = linspace(0, date35, Nt+1); dt = t(2)-t(1);

s35 = zeros(1,Nt); i35 = zeros(1,Nt); r35 = zeros(1,Nt);
s35(1) = population-i34(end); i35(1) = i34(end); r35(1) = r34(end);

b = 0.0000001250; bf35 = 0.001; k35 = 0.11; 
for j = 1:Nt
    s35(j+1) = s35(j) - b*bf35*s35(j)*i35(j)*dt;
    i35(j+1) = i35(j) + (b*bf35*s35(j)*i35(j) - (k35)*i35(j))*dt;
    r35(j+1) = r35(j) + k35*i35(j)*dt;
end
figure(35); clf;
hold on
plot(t, i35, 'b-', t, r35, 'g-');
plot(0:date35, MERS_I(35:36),'ro-', 0:date35, MERS_R(35:36),'mo-');
legend('Infected','Recovered',2);


%% day36~day37
clear Nt t dt b; format long;
date36 = datenum('25-june-2015') - datenum('24-june-2015');
Nt = sc*date36; t = linspace(0, date36, Nt+1); dt = t(2)-t(1);

s36 = zeros(1,Nt); i36 = zeros(1,Nt); r36 = zeros(1,Nt);
s36(1) = population-i35(end); i36(1) = i35(end); r36(1) = r35(end);

b = 0.0000001250; bf36 = 0.002; k36 = 0.12; 
for j = 1:Nt
    s36(j+1) = s36(j) - b*bf36*s36(j)*i36(j)*dt;
    i36(j+1) = i36(j) + (b*bf36*s36(j)*i36(j) - (k36)*i36(j))*dt;
    r36(j+1) = r36(j) + k36*i36(j)*dt;
end
figure(36); clf;
hold on
plot(t, i36, 'b-', t, r36, 'g-');
plot(0:date36, MERS_I(36:37),'ro-', 0:date36, MERS_R(36:37),'mo-');
legend('Infected','Recovered',2);


%% day37~day38
clear Nt t dt b; format long;
date37 = datenum('26-june-2015') - datenum('25-june-2015');
Nt = sc*date37; t = linspace(0, date37, Nt+1); dt = t(2)-t(1);

s37 = zeros(1,Nt); i37 = zeros(1,Nt); r37 = zeros(1,Nt);
s37(1) = population-i36(end); i37(1) = i36(end); r37(1) = r36(end);

b = 0.0000001250; bf37 = 0.003; k37 = 0.14; 
for j = 1:Nt
    s37(j+1) = s37(j) - b*bf37*s37(j)*i37(j)*dt;
    i37(j+1) = i37(j) + (b*bf37*s37(j)*i37(j) - (k37)*i37(j))*dt;
    r37(j+1) = r37(j) + k37*i37(j)*dt;
end
figure(37); clf;
hold on
plot(t, i37, 'b-', t, r37, 'g-');
plot(0:date37, MERS_I(37:38),'ro-', 0:date37, MERS_R(37:38),'mo-');
legend('Infected','Recovered',2);


%% day38~day39
clear Nt t dt b; format long;
date38 = datenum('27-june-2015') - datenum('26-june-2015');
Nt = sc*date38; t = linspace(0, date38, Nt+1); dt = t(2)-t(1);

s38 = zeros(1,Nt); i38 = zeros(1,Nt); r38 = zeros(1,Nt);
s38(1) = population-i37(end); i38(1) = i37(end); r38(1) = r37(end);

b = 0.0000001250; bf38 = 0.0003; k38 = 0.04; 
for j = 1:Nt
    s38(j+1) = s38(j) - b*bf38*s38(j)*i38(j)*dt;
    i38(j+1) = i38(j) + (b*bf38*s38(j)*i38(j) - (k38)*i38(j))*dt;
    r38(j+1) = r38(j) + k38*i38(j)*dt;
end
figure(38); clf;
hold on
plot(t, i38, 'b-', t, r38, 'g-');
plot(0:date38, MERS_I(38:39),'ro-', 0:date38, MERS_R(38:39),'mo-');
legend('Infected','Recovered',2);

%% day39~day40
clear Nt t dt b; format long;
date39 = datenum('28-june-2015') - datenum('27-june-2015');
Nt = sc*date39; t = linspace(0, date39, Nt+1); dt = t(2)-t(1);

s39 = zeros(1,Nt); i39 = zeros(1,Nt); r39 = zeros(1,Nt);
s39(1) = population-i38(end); i39(1) = i38(end); r39(1) = r38(end);

b = 0.0000001250; bf39 = 0.0; k39 = 0.03; 
for j = 1:Nt
    s39(j+1) = s39(j) - b*bf39*s39(j)*i39(j)*dt;
    i39(j+1) = i39(j) + (b*bf39*s39(j)*i39(j) - (k39)*i39(j))*dt;
    r39(j+1) = r39(j) + k39*i39(j)*dt;
end
figure(39); clf;
hold on
plot(t, i39, 'b-', t, r39, 'g-');
plot(0:date39, MERS_I(39:40),'ro-', 0:date39, MERS_R(39:40),'mo-');
legend('Infected','Recovered',2);

%% day40~day41
clear Nt t dt b; format long;
date40 = datenum('29-june-2015') - datenum('28-june-2015');
Nt = sc*date40; t = linspace(0, date40, Nt+1); dt = t(2)-t(1);

s40 = zeros(1,Nt); i40 = zeros(1,Nt); r40 = zeros(1,Nt);
s40(1) = population-i39(end); i40(1) = i39(end); r40(1) = r39(end);

b = 0.0000001250; bf40 = 0.0; k40 = 0.05; 
for j = 1:Nt
    s40(j+1) = s40(j) - b*bf40*s40(j)*i40(j)*dt;
    i40(j+1) = i40(j) + (b*bf40*s40(j)*i40(j) - (k40)*i40(j))*dt;
    r40(j+1) = r40(j) + k40*i40(j)*dt;
end
figure(40); clf;
hold on
plot(t, i40, 'b-', t, r40, 'g-');
plot(0:date40, MERS_I(40:41),'ro-', 0:date40, MERS_R(40:41),'mo-');
legend('Infected','Recovered',2);

%% day41~day42
clear Nt t dt b; format long;
date41 = datenum('30-june-2015') - datenum('29-june-2015');
Nt = sc*date41; t = linspace(0, date41, Nt+1); dt = t(2)-t(1);

s41 = zeros(1,Nt); i41 = zeros(1,Nt); r41 = zeros(1,Nt);
s41(1) = population-i40(end); i41(1) = i40(end); r41(1) = r40(end);

b = 0.0000001250; bf41 = 0.0; k41 = 0.04; 
for j = 1:Nt
    s41(j+1) = s41(j) - b*bf41*s41(j)*i41(j)*dt;
    i41(j+1) = i41(j) + (b*bf41*s41(j)*i41(j) - (k41)*i41(j))*dt;
    r41(j+1) = r41(j) + k41*i41(j)*dt;
end
figure(41); clf;
hold on
plot(t, i41, 'b-', t, r41, 'g-');
plot(0:date41, MERS_I(41:42),'ro-', 0:date41, MERS_R(41:42),'mo-');
legend('Infected','Recovered',2);


%% day42~day43
clear Nt t dt b; format long;
date42 = datenum('1-july-2015') - datenum('30-june-2015');
Nt = sc*date42; t = linspace(0, date42, Nt+1); dt = t(2)-t(1);

s42 = zeros(1,Nt); i42 = zeros(1,Nt); r42 = zeros(1,Nt);
s42(1) = population-i41(end); i42(1) = i41(end); r42(1) = r41(end);

b = 0.0000001250; bf42 = 0.003; k42 = 0.1; 
for j = 1:Nt
    s42(j+1) = s42(j) - b*bf42*s42(j)*i42(j)*dt;
    i42(j+1) = i42(j) + (b*bf42*s42(j)*i42(j) - (k42)*i42(j))*dt;
    r42(j+1) = r42(j) + k42*i42(j)*dt;
end
figure(42); clf;
hold on
plot(t, i42, 'b-', t, r42, 'g-');
plot(0:date42, MERS_I(42:43),'ro-', 0:date42, MERS_R(42:43),'mo-');
legend('Infected','Recovered',2);

%% day43~day44
clear Nt t dt b; format long;
date43 = datenum('2-july-2015') - datenum('1-july-2015');
Nt = sc*date43; t = linspace(0, date43, Nt+1); dt = t(2)-t(1);

s43 = zeros(1,Nt); i43 = zeros(1,Nt); r43 = zeros(1,Nt);
s43(1) = population-i42(end); i43(1) = i42(end); r43(1) = r42(end);

b = 0.0000001250; bf43 = 0.0; k43 = 0.15; 
for j = 1:Nt
    s43(j+1) = s43(j) - b*bf43*s43(j)*i43(j)*dt;
    i43(j+1) = i43(j) + (b*bf43*s43(j)*i43(j) - (k43)*i43(j))*dt;
    r43(j+1) = r43(j) + k43*i43(j)*dt;
end
figure(43); clf;
hold on
plot(t, i43, 'b-', t, r43, 'g-');
plot(0:date43, MERS_I(43:44),'ro-', 0:date43, MERS_R(43:44),'mo-');
legend('Infected','Recovered',2);

%% day44~day45
clear Nt t dt b; format long;
date44 = datenum('3-july-2015') - datenum('2-july-2015');
Nt = sc*date44; t = linspace(0, date44, Nt+1); dt = t(2)-t(1);

s44 = zeros(1,Nt); i44 = zeros(1,Nt); r44 = zeros(1,Nt);
s44(1) = population-i43(end); i44(1) = i43(end); r44(1) = r43(end);

b = 0.0000001250; bf44 = 0.005; k44 = 0.05; 
for j = 1:Nt
    s44(j+1) = s44(j) - b*bf44*s44(j)*i44(j)*dt;
    i44(j+1) = i44(j) + (b*bf44*s44(j)*i44(j) - (k44)*i44(j))*dt;
    r44(j+1) = r44(j) + k44*i44(j)*dt;
end
figure(44); clf;
hold on
plot(t, i44, 'b-', t, r44, 'g-');
plot(0:date44, MERS_I(44:45),'ro-', 0:date44, MERS_R(44:45),'mo-');
legend('Infected','Recovered',2);


%% day45~day46
clear Nt t dt b; format long;
date45 = datenum('4-july-2015') - datenum('3-july-2015');
Nt = sc*date45; t = linspace(0, date45, Nt+1); dt = t(2)-t(1);

s45 = zeros(1,Nt); i45 = zeros(1,Nt); r45 = zeros(1,Nt);
s45(1) = population-i44(end); i45(1) = i44(end); r45(1) = r44(end);

b = 0.0000001250; bf45 = 0.008; k45 = 0.15; 
for j = 1:Nt
    s45(j+1) = s45(j) - b*bf45*s45(j)*i45(j)*dt;
    i45(j+1) = i45(j) + (b*bf45*s45(j)*i45(j) - (k45)*i45(j))*dt;
    r45(j+1) = r45(j) + k45*i45(j)*dt;
end
figure(45); clf;
hold on
plot(t, i45, 'b-', t, r45, 'g-');
plot(0:date45, MERS_I(45:46),'ro-', 0:date45, MERS_R(45:46),'mo-');
legend('Infected','Recovered',2);

%% day46~day47
clear Nt t dt b; format long;
date46 = datenum('5-july-2015') - datenum('4-july-2015');
Nt = sc*date46; t = linspace(0, date46, Nt+1); dt = t(2)-t(1);

s46 = zeros(1,Nt); i46 = zeros(1,Nt); r46 = zeros(1,Nt);
s46(1) = population-i45(end); i46(1) = i45(end); r46(1) = r45(end);

b = 0.0000001250; bf46 = 0.0; k46 = 0.02; 
for j = 1:Nt
    s46(j+1) = s46(j) - b*bf46*s46(j)*i46(j)*dt;
    i46(j+1) = i46(j) + (b*bf46*s46(j)*i46(j) - (k46)*i46(j))*dt;
    r46(j+1) = r46(j) + k46*i46(j)*dt;
end
figure(46); clf;
hold on
plot(t, i46, 'b-', t, r46, 'g-');
plot(0:date46, MERS_I(46:47),'ro-', 0:date46, MERS_R(46:47),'mo-');
legend('Infected','Recovered',2);


%% day47~day48
clear Nt t dt b; format long;
date47 = datenum('6-july-2015') - datenum('5-july-2015');
Nt = sc*date47; t = linspace(0, date47, Nt+1); dt = t(2)-t(1);

s47 = zeros(1,Nt); i47 = zeros(1,Nt); r47 = zeros(1,Nt);
s47(1) = population-i46(end); i47(1) = i46(end); r47(1) = r46(end);

b = 0.0000001250; bf47 = 0.0; k47 = 0.02; 
for j = 1:Nt
    s47(j+1) = s47(j) - b*bf47*s47(j)*i47(j)*dt;
    i47(j+1) = i47(j) + (b*bf47*s47(j)*i47(j) - (k47)*i47(j))*dt;
    r47(j+1) = r47(j) + k47*i47(j)*dt;
end
figure(47); clf;
hold on
plot(t, i47, 'b-', t, r47, 'g-');
plot(0:date47, MERS_I(47:48),'ro-', 0:date47, MERS_R(47:48),'mo-');
legend('Infected','Recovered',2);



%% day48~day49
clear Nt t dt b; format long;
date48 = datenum('7-july-2015') - datenum('6-july-2015');
Nt = sc*date48; t = linspace(0, date48, Nt+1); dt = t(2)-t(1);

s48 = zeros(1,Nt); i48 = zeros(1,Nt); r48 = zeros(1,Nt);
s48(1) = population-i47(end); i48(1) = i47(end); r48(1) = r47(end);

b = 0.0000001250; bf48 = 0.0; k48 = 0.06; 
for j = 1:Nt
    s48(j+1) = s48(j) - b*bf48*s48(j)*i48(j)*dt;
    i48(j+1) = i48(j) + (b*bf48*s48(j)*i48(j) - (k48)*i48(j))*dt;
    r48(j+1) = r48(j) + k48*i48(j)*dt;
end
figure(48); clf;
hold on
plot(t, i48, 'b-', t, r48, 'g-');
plot(0:date48, MERS_I(48:49),'ro-', 0:date48, MERS_R(48:49),'mo-');
legend('Infected','Recovered',2);

%% day49~day50
clear Nt t dt b; format long;
date49 = datenum('8-july-2015') - datenum('7-july-2015');
Nt = sc*date49; t = linspace(0, date49, Nt+1); dt = t(2)-t(1);

s49 = zeros(1,Nt); i49 = zeros(1,Nt); r49 = zeros(1,Nt);
s49(1) = population-i48(end); i49(1) = i48(end); r49(1) = r48(end);

b = 0.0000001250; bf49 = 0.0; k49 = 0.06; 
for j = 1:Nt
    s49(j+1) = s49(j) - b*bf49*s49(j)*i49(j)*dt;
    i49(j+1) = i49(j) + (b*bf49*s49(j)*i49(j) - (k49)*i49(j))*dt;
    r49(j+1) = r49(j) + k49*i49(j)*dt;
end
figure(49); clf;
hold on
plot(t, i49, 'b-', t, r49, 'g-');
plot(0:date49, MERS_I(49:50),'ro-', 0:date49, MERS_R(49:50),'mo-');
legend('Infected','Recovered',2);

%% day50~day51
clear Nt t dt b; format long;
date50 = datenum('9-july-2015') - datenum('8-july-2015');
Nt = sc*date50; t = linspace(0, date50, Nt+1); dt = t(2)-t(1);

s50 = zeros(1,Nt); i50 = zeros(1,Nt); r50 = zeros(1,Nt);
s50(1) = population-i49(end); i50(1) = i49(end); r50(1) = r49(end);

b = 0.0000001250; bf50 = 0.0; k50 = 0.2; 
for j = 1:Nt
    s50(j+1) = s50(j) - b*bf50*s50(j)*i50(j)*dt;
    i50(j+1) = i50(j) + (b*bf50*s50(j)*i50(j) - (k50)*i50(j))*dt;
    r50(j+1) = r50(j) + k50*i50(j)*dt;
end
figure(50); clf;
hold on
plot(t, i50, 'b-', t, r50, 'g-');
plot(0:date50, MERS_I(50:51),'ro-', 0:date50, MERS_R(50:51),'mo-');
legend('Infected','Recovered',2);


%% day51~day52
clear Nt t dt b; format long;
date51 = datenum('10-july-2015') - datenum('9-july-2015');
Nt = sc*date51; t = linspace(0, date51, Nt+1); dt = t(2)-t(1);

s51 = zeros(1,Nt); i51 = zeros(1,Nt); r51 = zeros(1,Nt);
s51(1) = population-i50(end); i51(1) = i50(end); r51(1) = r50(end);

b = 0.0000001250; bf51 = 0.0; k51 = 0.16; 
for j = 1:Nt
    s51(j+1) = s51(j) - b*bf51*s51(j)*i51(j)*dt;
    i51(j+1) = i51(j) + (b*bf51*s51(j)*i51(j) - (k51)*i51(j))*dt;
    r51(j+1) = r51(j) + k51*i51(j)*dt;
end
figure(51); clf;
hold on
plot(t, i51, 'b-', t, r51, 'g-');
plot(0:date51, MERS_I(51:52),'ro-', 0:date51, MERS_R(51:52),'mo-');
legend('Infected','Recovered',2);

%% day52~day53
clear Nt t dt b; format long;
date52 = datenum('11-july-2015') - datenum('10-july-2015');
Nt = sc*date52; t = linspace(0, date52, Nt+1); dt = t(2)-t(1);

s52 = zeros(1,Nt); i52 = zeros(1,Nt); r52 = zeros(1,Nt);
s52(1) = population-i51(end); i52(1) = i51(end); r52(1) = r51(end);

b = 0.0000001250; bf52 = 0.0; k52 = 0.07; 
for j = 1:Nt
    s52(j+1) = s52(j) - b*bf52*s52(j)*i52(j)*dt;
    i52(j+1) = i52(j) + (b*bf52*s52(j)*i52(j) - (k52)*i52(j))*dt;
    r52(j+1) = r52(j) + k52*i52(j)*dt;
end
figure(52); clf;
hold on
plot(t, i52, 'b-', t, r52, 'g-');
plot(0:date52, MERS_I(52:53),'ro-', 0:date52, MERS_R(52:53),'mo-');
legend('Infected','Recovered',2);

%% day53~day54
clear Nt t dt b; format long;
date53 = datenum('12-july-2015') - datenum('11-july-2015');
Nt = sc*date53; t = linspace(0, date53, Nt+1); dt = t(2)-t(1);

s53 = zeros(1,Nt); i53 = zeros(1,Nt); r53 = zeros(1,Nt);
s53(1) = population-i52(end); i53(1) = i52(end); r53(1) = r52(end);

b = 0.0000001250; bf53 = 0.0; k53 = 0.02; 
for j = 1:Nt
    s53(j+1) = s53(j) - b*bf53*s53(j)*i53(j)*dt;
    i53(j+1) = i53(j) + (b*bf53*s53(j)*i53(j) - (k53)*i53(j))*dt;
    r53(j+1) = r53(j) + k53*i53(j)*dt;
end
figure(53); clf;
hold on
plot(t, i53, 'b-', t, r53, 'g-');
plot(0:date53, MERS_I(53:54),'ro-', 0:date53, MERS_R(53:54),'mo-');
legend('Infected','Recovered',2);


%% day54~day55
clear Nt t dt b; format long;
date54 = datenum('13-july-2015') - datenum('12-july-2015');
Nt = sc*date54; t = linspace(0, date54, Nt+1); dt = t(2)-t(1);

s54 = zeros(1,Nt); i54 = zeros(1,Nt); r54 = zeros(1,Nt);
s54(1) = population-i53(end); i54(1) = i53(end); r54(1) = r53(end);

b = 0.0000001250; bf54 = 0.0; k54 = 0.04; 
for j = 1:Nt
    s54(j+1) = s54(j) - b*bf54*s54(j)*i54(j)*dt;
    i54(j+1) = i54(j) + (b*bf54*s54(j)*i54(j) - (k54)*i54(j))*dt;
    r54(j+1) = r54(j) + k54*i54(j)*dt;
end
figure(54); clf;
hold on
plot(t, i54, 'b-', t, r54, 'g-');
plot(0:date54, MERS_I(54:55),'ro-', 0:date54, MERS_R(54:55),'mo-');
legend('Infected','Recovered',2);



%% day55~day56
clear Nt t dt b; format long;
date55 = datenum('14-july-2015') - datenum('13-july-2015');
Nt = sc*date55; t = linspace(0, date55, Nt+1); dt = t(2)-t(1);

s55 = zeros(1,Nt); i55 = zeros(1,Nt); r55 = zeros(1,Nt);
s55(1) = population-i54(end); i55(1) = i54(end); r55(1) = r54(end);

b = 0.0000001250; bf55 = 0.0; k55 = 0.05; 
for j = 1:Nt
    s55(j+1) = s55(j) - b*bf55*s55(j)*i55(j)*dt;
    i55(j+1) = i55(j) + (b*bf55*s55(j)*i55(j) - (k55)*i55(j))*dt;
    r55(j+1) = r55(j) + k55*i55(j)*dt;
end
figure(55); clf;
hold on
plot(t, i55, 'b-', t, r55, 'g-');
plot(0:date55, MERS_I(55:56),'ro-', 0:date55, MERS_R(55:56),'mo-');
legend('Infected','Recovered',2);

%% day56~day57
clear Nt t dt b; format long;
date56 = datenum('15-july-2015') - datenum('14-july-2015');
Nt = sc*date56; t = linspace(0, date56, Nt+1); dt = t(2)-t(1);

s56 = zeros(1,Nt); i56 = zeros(1,Nt); r56 = zeros(1,Nt);
s56(1) = population-i55(end); i56(1) = i55(end); r56(1) = r55(end);

b = 0.0000001250; bf56 = 0.0; k56 = 0.08; 
for j = 1:Nt
    s56(j+1) = s56(j) - b*bf56*s56(j)*i56(j)*dt;
    i56(j+1) = i56(j) + (b*bf56*s56(j)*i56(j) - (k56)*i56(j))*dt;
    r56(j+1) = r56(j) + k56*i56(j)*dt;
end
figure(56); clf;
hold on
plot(t, i56, 'b-', t, r56, 'g-');
plot(0:date56, MERS_I(56:57),'ro-', 0:date56, MERS_R(56:57),'mo-');
legend('Infected','Recovered',2);


%% All days
figure(100); clf; LW = 2.5; t = linspace(0, 1, 101);
plot(t, i1, 'b-', t, r1, 'g-','linewidth',LW);
hold on
plot(t+1, i2, 'b-', t+1, r2, 'g-','linewidth',LW);
plot(t+2, i3, 'b-', t+2, r3, 'g-','linewidth',LW);
plot(t+3, i4, 'b-', t+3, r4, 'g-','linewidth',LW);
plot(t+4, i5, 'b-', t+4, r5, 'g-','linewidth',LW);
plot(t+5, i6, 'b-', t+5, r6, 'g-','linewidth',LW);
plot(t+6, i7, 'b-', t+6, r7, 'g-','linewidth',LW);
plot(t+7, i8, 'b-', t+7, r8, 'g-','linewidth',LW);
plot(t+8, i9, 'b-', t+8, r9, 'g-','linewidth',LW);
plot(t+9, i10, 'b-', t+9, r10, 'g-','linewidth',LW);
plot(t+10, i11, 'b-', t+10, r11, 'g-','linewidth',LW);
plot(t+11, i12, 'b-', t+11, r12, 'g-','linewidth',LW);
plot(t+12, i13, 'b-', t+12, r13, 'g-','linewidth',LW);
plot(t+13, i14, 'b-', t+13, r14, 'g-','linewidth',LW);
plot(t+14, i15, 'b-', t+14, r15, 'g-','linewidth',LW);
plot(t+15, i16, 'b-', t+15, r16, 'g-','linewidth',LW);
plot(t+16, i17, 'b-', t+16, r17, 'g-','linewidth',LW);
plot(t+17, i18, 'b-', t+17, r18, 'g-','linewidth',LW);
plot(t+18, i19, 'b-', t+18, r19, 'g-','linewidth',LW);
plot(t+19, i20, 'b-', t+19, r20, 'g-','linewidth',LW);
plot(t+20, i21, 'b-', t+20, r21, 'g-','linewidth',LW);
plot(t+21, i22, 'b-', t+21, r22, 'g-','linewidth',LW);
plot(t+22, i23, 'b-', t+22, r23, 'g-','linewidth',LW);
plot(t+23, i24, 'b-', t+23, r24, 'g-','linewidth',LW);
plot(t+24, i25, 'b-', t+24, r25, 'g-','linewidth',LW);
plot(t+25, i26, 'b-', t+25, r26, 'g-','linewidth',LW);
plot(t+26, i27, 'b-', t+26, r27, 'g-','linewidth',LW);
plot(t+27, i28, 'b-', t+27, r28, 'g-','linewidth',LW);
plot(t+28, i29, 'b-', t+28, r29, 'g-','linewidth',LW);
plot(t+29, i30, 'b-', t+29, r30, 'g-','linewidth',LW);
plot(t+30, i31, 'b-', t+30, r31, 'g-','linewidth',LW);
plot(t+31, i32, 'b-', t+31, r32, 'g-','linewidth',LW);
plot(t+32, i33, 'b-', t+32, r33, 'g-','linewidth',LW);
plot(t+33, i34, 'b-', t+33, r34, 'g-','linewidth',LW);
plot(t+34, i35, 'b-', t+34, r35, 'g-','linewidth',LW);
plot(t+35, i36, 'b-', t+35, r36, 'g-','linewidth',LW);
plot(t+36, i37, 'b-', t+36, r37, 'g-','linewidth',LW);
plot(t+37, i38, 'b-', t+37, r38, 'g-','linewidth',LW);
plot(t+38, i39, 'b-', t+38, r39, 'g-','linewidth',LW);
plot(t+39, i40, 'b-', t+39, r40, 'g-','linewidth',LW);
plot(t+40, i41, 'b-', t+40, r41, 'g-','linewidth',LW);
plot(t+41, i42, 'b-', t+41, r42, 'g-','linewidth',LW);
plot(t+42, i43, 'b-', t+42, r43, 'g-','linewidth',LW);
plot(t+43, i44, 'b-', t+43, r44, 'g-','linewidth',LW);
plot(t+44, i45, 'b-', t+44, r45, 'g-','linewidth',LW);
plot(t+45, i46, 'b-', t+45, r46, 'g-','linewidth',LW);
plot(t+46, i47, 'b-', t+46, r47, 'g-','linewidth',LW);
plot(t+47, i48, 'b-', t+47, r48, 'g-','linewidth',LW);
plot(t+48, i49, 'b-', t+48, r49, 'g-','linewidth',LW);
plot(t+49, i50, 'b-', t+49, r50, 'g-','linewidth',LW);
plot(t+50, i51, 'b-', t+50, r51, 'g-','linewidth',LW);
plot(t+51, i52, 'b-', t+51, r52, 'g-','linewidth',LW);
plot(t+52, i53, 'b-', t+52, r53, 'g-','linewidth',LW);
plot(t+53, i54, 'b-', t+53, r54, 'g-','linewidth',LW);
plot(t+54, i55, 'b-', t+54, r55, 'g-','linewidth',LW);
plot(t+55, i56, 'b-', t+55, r56, 'g-','linewidth',LW);
plot(0:date, MERS_I,'ro-', 0:date, MERS_R,'mo-');
legend('Infected','Recovered' ,2);


%% b value
figure(200); clf
plot(t, b*bf1, 'b-');
hold on
plot(t+1, b*bf2, 'b-');
plot(t+2, b*bf3, 'b-');
plot(t+3, b*bf4, 'b-');
plot(t+4, b*bf5, 'b-');
plot(t+5, b*bf6, 'b-');
plot(t+6, b*bf7, 'b-');
plot(t+7, b*bf8, 'b-');
plot(t+8, b*bf9, 'b-');
plot(t+9, b*bf10, 'b-');
plot(t+10, b*bf11, 'b-');
plot(t+11, b*bf12, 'b-');
plot(t+12, b*bf13, 'b-');
plot(t+13, b*bf14, 'b-');
plot(t+14, b*bf15, 'b-');
plot(t+15, b*bf16, 'b-');
plot(t+16, b*bf17, 'b-');
plot(t+17, b*bf18, 'b-');
plot(t+18, b*bf19, 'b-');
plot(t+19, b*bf20, 'b-');
plot(t+20, b*bf21, 'b-');
plot(t+21, b*bf22, 'b-');
plot(t+22, b*bf23, 'b-');
plot(t+23, b*bf24, 'b-');
plot(t+24, b*bf25, 'b-');
plot(t+25, b*bf26, 'b-');
plot(t+26, b*bf27, 'b-');
plot(t+27, b*bf28, 'b-');
plot(t+28, b*bf29, 'b-');
plot(t+29, b*bf30, 'b-');
plot(t+30, b*bf31, 'b-');
plot(t+31, b*bf32, 'b-');
plot(t+32, b*bf33, 'b-');
plot(t+33, b*bf34, 'b-');
plot(t+34, b*bf35, 'b-');
plot(t+35, b*bf36, 'b-');
plot(t+36, b*bf37, 'b-');
plot(t+37, b*bf38, 'b-');
plot(t+38, b*bf39, 'b-');
plot(t+39, b*bf40, 'b-');
plot(t+40, b*bf41, 'b-');
plot(t+41, b*bf42, 'b-');
plot(t+42, b*bf43, 'b-');
plot(t+43, b*bf44, 'b-');
plot(t+44, b*bf45, 'b-');
plot(t+45, b*bf46, 'b-');
plot(t+46, b*bf47, 'b-');
plot(t+47, b*bf48, 'b-');
plot(t+48, b*bf49, 'b-');
plot(t+49, b*bf50, 'b-');
plot(t+50, b*bf51, 'b-');
plot(t+51, b*bf52, 'b-');
plot(t+52, b*bf53, 'b-');
plot(t+53, b*bf54, 'b-');
plot(t+54, b*bf55, 'b-');
plot(t+55, b*bf56, 'b-');
legend('beta')

%% k value
figure(300); clf
plot(t, k1, 'b-');
hold on
plot(t+1, k2, 'b-');
plot(t+2, k3, 'b-');
plot(t+3, k4, 'b-');
plot(t+4, k5, 'b-');
plot(t+5, k6, 'b-');
plot(t+6, k7, 'b-');
plot(t+7, k8, 'b-');
plot(t+8, k9, 'b-');
plot(t+9, k10, 'b-');
plot(t+10, k11, 'b-');
plot(t+11, k12, 'b-');
plot(t+12, k13, 'b-');
plot(t+13, k14, 'b-');
plot(t+14, k15, 'b-');
plot(t+15, k16, 'b-');
plot(t+16, k17, 'b-');
plot(t+17, k18, 'b-');
plot(t+18, k19, 'b-');
plot(t+19, k20, 'b-');
plot(t+20, k21, 'b-');
plot(t+21, k22, 'b-');
plot(t+22, k23, 'b-');
plot(t+23, k24, 'b-');
plot(t+24, k25, 'b-');
plot(t+25, k26, 'b-');
plot(t+26, k27, 'b-');
plot(t+27, k28, 'b-');
plot(t+28, k29, 'b-');
plot(t+29, k30, 'b-');
plot(t+30, k31, 'b-');
plot(t+31, k32, 'b-');
plot(t+32, k33, 'b-');
plot(t+33, k34, 'b-');
plot(t+34, k35, 'b-');
plot(t+35, k36, 'b-');
plot(t+36, k37, 'b-');
plot(t+37, k38, 'b-');
plot(t+38, k39, 'b-');
plot(t+39, k40, 'b-');
plot(t+40, k41, 'b-');
plot(t+41, k42, 'b-');
plot(t+42, k43, 'b-');
plot(t+43, k44, 'b-');
plot(t+44, k45, 'b-');
plot(t+45, k46, 'b-');
plot(t+46, k47, 'b-');
plot(t+47, k48, 'b-');
plot(t+48, k49, 'b-');
plot(t+49, k50, 'b-');
plot(t+50, k51, 'b-');
plot(t+51, k52, 'b-');
plot(t+52, k53, 'b-');
plot(t+53, k54, 'b-');
plot(t+54, k55, 'b-');
plot(t+55, k56, 'b-');
legend('gamma')

%% R0
P = population;
R0_1 = (b*bf1*P)/k1;
R0_2 = (b*bf2*P)/k2;
R0_3 = (b*bf3*P)/k3;
R0_4 = (b*bf4*P)/k4;
R0_5 = (b*bf5*P)/k5;
R0_6 = (b*bf6*P)/k6;
R0_7 = (b*bf7*P)/k7;
R0_8 = (b*bf8*P)/k8;
R0_9 = (b*bf9*P)/k9;
R0_10 = (b*bf10*P)/k10;
R0_11 = (b*bf11*P)/k11;
R0_12 = (b*bf12*P)/k12;
R0_13 = (b*bf13*P)/k13;
R0_14 = (b*bf14*P)/k14;
R0_15 = (b*bf15*P)/k15;
R0_16 = (b*bf16*P)/k16;
R0_17 = (b*bf17*P)/k17;
R0_18 = (b*bf18*P)/k18;
R0_19 = (b*bf19*P)/k19;
R0_20 = (b*bf20*P)/k20;
R0_21 = (b*bf21*P)/k21;
R0_22 = (b*bf22*P)/k22;
R0_23 = (b*bf23*P)/k23;
R0_24 = (b*bf24*P)/k24;
R0_25 = (b*bf25*P)/k25;
R0_26 = (b*bf26*P)/k26;
R0_27 = (b*bf27*P)/k27;
R0_28 = (b*bf28*P)/k28;
R0_29 = (b*bf29*P)/k29;
R0_30 = (b*bf30*P)/k30;
R0_31 = (b*bf31*P)/k31;
R0_32 = (b*bf32*P)/k32;
R0_33 = (b*bf33*P)/k33;
R0_34 = (b*bf34*P)/k34;
R0_35 = (b*bf35*P)/k35;
R0_36 = (b*bf36*P)/k36;
R0_37 = (b*bf37*P)/k37;
R0_38 = (b*bf38*P)/k38;
R0_39 = (b*bf39*P)/k39;
R0_40 = (b*bf40*P)/k40;
R0_41 = (b*bf41*P)/k41;
R0_42 = (b*bf42*P)/k42;
R0_43 = (b*bf43*P)/k43;
R0_44 = (b*bf44*P)/k44;
R0_45 = (b*bf45*P)/k45;
R0_46 = (b*bf46*P)/k46;
R0_47 = (b*bf47*P)/k47;
R0_48 = (b*bf48*P)/k48;
R0_49 = (b*bf49*P)/k49;
R0_50 = (b*bf50*P)/k50;
R0_51 = (b*bf51*P)/k51;
R0_52 = (b*bf52*P)/k52;
R0_53 = (b*bf53*P)/k53;
R0_54 = (b*bf54*P)/k54;
R0_55 = (b*bf55*P)/k55;
R0_56 = (b*bf56*P)/k56;


figure(400); clf;
plot(t,R0_1,'b-');hold on
plot(t+1,R0_2,'b-');
plot(t+2,R0_3,'b-');
plot(t+3,R0_4,'b-');
plot(t+4,R0_5,'b-');
plot(t+5,R0_6,'b-');
plot(t+6,R0_7,'b-');
plot(t+7,R0_8,'b-');
plot(t+8,R0_9,'b-');
plot(t+9,R0_10,'b-');
plot(t+10,R0_11,'b-');
plot(t+11,R0_12,'b-');
plot(t+12,R0_13,'b-');
plot(t+13,R0_14,'b-');
plot(t+14,R0_15,'b-');
plot(t+15,R0_16,'b-');
plot(t+16,R0_17,'b-');
plot(t+17,R0_18,'b-');
plot(t+18,R0_19,'b-');
plot(t+19,R0_20,'b-');
plot(t+20,R0_21,'b-');
plot(t+21,R0_22,'b-');
plot(t+22,R0_23,'b-');
plot(t+23,R0_24,'b-');
plot(t+24,R0_25,'b-');
plot(t+25,R0_26,'b-');
plot(t+26,R0_27,'b-');
plot(t+27,R0_28,'b-');
plot(t+28,R0_29,'b-');
plot(t+29,R0_30,'b-');
plot(t+30,R0_31,'b-');
plot(t+31,R0_32,'b-');
plot(t+32,R0_33,'b-');
plot(t+33,R0_34,'b-');
plot(t+34,R0_35,'b-');
plot(t+35,R0_36,'b-');
plot(t+36,R0_37,'b-');
plot(t+37,R0_38,'b-');
plot(t+38,R0_39,'b-');
plot(t+39,R0_40,'b-');
plot(t+40,R0_41,'b-');
plot(t+41,R0_42,'b-');
plot(t+42,R0_43,'b-');
plot(t+43,R0_44,'b-');
plot(t+44,R0_45,'b-');
plot(t+45,R0_46,'b-');
plot(t+46,R0_47,'b-');
plot(t+47,R0_48,'b-');
plot(t+48,R0_49,'b-');
plot(t+49,R0_50,'b-');
plot(t+50,R0_51,'b-');
plot(t+51,R0_52,'b-');
plot(t+52,R0_53,'b-');
plot(t+53,R0_54,'b-');
plot(t+54,R0_55,'b-');
plot(t+55,R0_56,'b-');
legend('R0',1)


A=[R0_1 R0_2 R0_3 R0_4 R0_5 R0_6 R0_7 R0_8 R0_9 R0_10 R0_11 R0_12 R0_13 R0_14 R0_15 R0_16 R0_17 R0_18 R0_19 R0_20 R0_21 R0_22 R0_23 R0_24 R0_25 R0_26 R0_27 R0_28 R0_29 R0_30 R0_31 R0_32 R0_33 R0_34 R0_35 R0_36 R0_37 R0_38 R0_39 R0_40 R0_41 R0_42 R0_43 R0_44 R0_45 R0_46 R0_47 R0_48 R0_49 R0_50 R0_51 R0_52 R0_53 R0_54 R0_55 R0_56]';

figure(10000)
for i=1:56
    plot(i+(0:0.01:1)-1,A(i),'r-')
    hold on
end


%% important graphs

A=[R0_1 R0_2 R0_3 R0_4 R0_5 R0_6 R0_7 R0_8 R0_9 R0_10 R0_11 R0_12 R0_13 R0_14 R0_15 R0_16 R0_17 R0_18 R0_19 R0_20 R0_21 R0_22 R0_23 R0_24 R0_25 R0_26 R0_27 R0_28 R0_29 R0_30 R0_31 R0_32 R0_33 R0_34 R0_35 R0_36 R0_37 R0_38 R0_39 R0_40 R0_41 R0_42 R0_43 R0_44 R0_45 R0_46 R0_47 R0_48 R0_49 R0_50 R0_51 R0_52 R0_53 R0_54 R0_55 R0_56]';
B=[b*bf1 b*bf2 b*bf3 b*bf4 b*bf5 b*bf6 b*bf7 b*bf8 b*bf9 b*bf10 b*bf11 b*bf12 b*bf13 b*bf14 b*bf15 b*bf16 b*bf17 b*bf18 b*bf19 b*bf20 b*bf21 b*bf22 b*bf23 b*bf24 b*bf25 b*bf26 b*bf27 b*bf28 b*bf29 b*bf30 b*bf31 b*bf32 b*bf33 b*bf34 b*bf35 b*bf36 b*bf37 b*bf38 b*bf39 b*bf40 b*bf41 b*bf42 b*bf43 b*bf44 b*bf45 b*bf46 b*bf47 b*bf48 b*bf49 b*bf50 b*bf51 b*bf52 b*bf53 b*bf54 b*bf55 b*bf56]
C=[k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18 k19 k20 k21 k22 k23 k24 k25 k26 k27 k28 k29 k30 k31 k32 k33 k34 k35 k36 k37 k38 k39 k40 k41 k42 k43 k44 k45 k46 k47 k48 k49 k50 k51 k52 k53 k54 k55 k56]



figure(1000)
plot(0:55,A,'k-');hold on
plot(0:55,1,'k-')
legend('R0',2)


figure(2000)
plot(0:55,B,'b-')
legend('beta',2)

figure(3000)
plot(0:55,C,'r-')
legend('gamma',2)