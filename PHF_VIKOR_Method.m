clear;
clc;
Original_Matrix=...
    {{0.2,[0.6,0.8]},{[0.6,0.8],[0.1,0.3]},{[0.4,0.5],[0.2,0.3]},{[0.5,0.6],0.5};...
    {[0.3,0.5],0.5},{[0.3,0.4,0.6],0.6},{0.2,[0.5,0.6]},{[0.5,0.7],[0.4,0.5]};...
    {0.7,[0.3,0.4,0.5]},{[0.5,0.6,0.7],0.7},{[0.7,0.8],[0.4,0.5]},{[0.3,0.4],0.4};...
    {[0.4,0.5],[0.2,0.4]},{0.6,[0.2,0.3]},{[0.3,0.5,0.6],0.5},{[0.7,0.8],[0.2,0.3]}}
Standard_Matrix=Original_Matrix;
Standard_Matrix_Size=size(Standard_Matrix);
Standard_Matrix_Row=Standard_Matrix_Size(1);
Standard_Matrix_Column=Standard_Matrix_Size(2);
for i=[1:Standard_Matrix_Row]
    Exchange_Unit=Standard_Matrix{i,1}(1);
    Standard_Matrix{i,1}(1)=Standard_Matrix{i,1}(2);
    Standard_Matrix{i,1}(2)=Exchange_Unit;
end
display(Standard_Matrix)
V_Membership_Token=[ ];
for i=[1:Standard_Matrix_Row]
    for j=[1:Standard_Matrix_Column]
        if max(find(Standard_Matrix{i,j}{1}))==1
            V_Membership_Token(end+1)=Standard_Matrix{i,j}{1}(1);
        elseif  max(find(Standard_Matrix{i,j}{1}))==2
            V_Membership_Token(end+1)=(1/2)*sqrt((Standard_Matrix{i,j}{1}(2)-Standard_Matrix{i,j}{1}(1))^2);
        else if max(find(Standard_Matrix{i,j}{1}))==3
                V_Membership_Token(end+1)=(1/3)*sqrt((Standard_Matrix{i,j}{1}(3)-Standard_Matrix{i,j}{1}(2))^2+...
                    (Standard_Matrix{i,j}{1}(3)-Standard_Matrix{i,j}{1}(1)^2+(Standard_Matrix{i,j}{1}(2)-Standard_Matrix{i,j}{1}(1))^2));
            end
        end
    end
end
V_Membership=[ ];
for i=[1:Standard_Matrix_Column:Standard_Matrix_Row*Standard_Matrix_Column]
    Row=[V_Membership_Token(i),V_Membership_Token(i+1),V_Membership_Token(i+2),V_Membership_Token(i+3)];
    V_Membership=[V_Membership;Row];
end
display(V_Membership)
Delta=0.1
% Delta is the coefficient of score.
G_Membership_Token=[ ];
for i=[1:Standard_Matrix_Row]
    for j=[1:Standard_Matrix_Column]
        G_Membership_Token(end+1)=(1/max(find(Standard_Matrix{i,j}{1})))*sum((Standard_Matrix{i,j}{1}).^0.1)*(1/0.1);
    end
end
G_Membership=[ ];
for i=[1:Standard_Matrix_Row:Standard_Matrix_Row*Standard_Matrix_Column]
    Row=[G_Membership_Token(i),G_Membership_Token(i+1),G_Membership_Token(i+2),G_Membership_Token(i+3)];
    G_Membership=[G_Membership;Row];
end
display(G_Membership)
S_Membership=V_Membership./(V_Membership+G_Membership)
E_Membership_Token=[ ];
for i=[1:Standard_Matrix_Column]
    for j=[1:Standard_Matrix_Column]
        if max(find(Standard_Matrix{i,j}{1}))==1
            E_Membership_Token(end+1)=Standard_Matrix{i,j}{1}(1);
        elseif  max(find(Standard_Matrix{i,j}{1}))==2
            E_Membership_Token(end+1)=(1/2)*abs(Standard_Matrix{i,j}{1}(2)-Standard_Matrix{i,j}{1}(1));
        else if max(find(Standard_Matrix{i,j}{1}))==3
                E_Membership_Token(end+1)=(1/3)*(abs(Standard_Matrix{i,j}{1}(1)-Standard_Matrix{i,j}{1}(2))+...
                    abs(Standard_Matrix{i,j}{1}(2)-Standard_Matrix{i,j}{1}(3)));
            end
        end
    end
end
E_Membership=[ ];
for i=[1:Standard_Matrix_Row:Standard_Matrix_Row*Standard_Matrix_Column]
    Row=[E_Membership_Token(i),E_Membership_Token(i+1),E_Membership_Token(i+2),E_Membership_Token(i+3)];
    E_Membership=[E_Membership;Row];
end
display(E_Membership)
Reverse_Standard_Matrix=Standard_Matrix;
for i=[1:Standard_Matrix_Row]
    for j=[1:Standard_Matrix_Column]
        for k=[1:max(find(Standard_Matrix{i,j}{2}))]
            Reverse_Standard_Matrix{i,j}{2}(k)=1-Standard_Matrix{i,j}{2}(k);
        end
    end
end
V_NonMembership_Token=[ ];
for i=[1:Standard_Matrix_Row]
    for j=[1:Standard_Matrix_Column]
        if max(find(Reverse_Standard_Matrix{i,j}{1}))==1
            V_NonMembership_Token(end+1)=Reverse_Standard_Matrix{i,j}{1}(1);
        elseif  max(find(Reverse_Standard_Matrix{i,j}{1}))==2
            V_NonMembership_Token(end+1)=(1/2)*sqrt((Reverse_Standard_Matrix{i,j}{1}(2)-Reverse_Standard_Matrix{i,j}{1}(1))^2);
        else if max(find(Reverse_Standard_Matrix{i,j}{1}))==3
                V_NonMembership_Token(end+1)=(1/3)*sqrt((Reverse_Standard_Matrix{i,j}{1}(3)-Reverse_Standard_Matrix{i,j}{1}(2))^2+...
                    (Reverse_Standard_Matrix{i,j}{1}(3)-Reverse_Standard_Matrix{i,j}{1}(1)^2+...
                    (Reverse_Standard_Matrix{i,j}{1}(2)-Reverse_Standard_Matrix{i,j}{1}(1))^2));
            end
        end
    end
end
V_NonMembership=[ ];
for i=[1:Standard_Matrix_Row:Standard_Matrix_Row*Standard_Matrix_Column]
    Row=[V_NonMembership_Token(i),V_NonMembership_Token(i+1),V_NonMembership_Token(i+2),V_NonMembership_Token(i+3)];
    V_NonMembership=[V_NonMembership;Row];
end
display(V_NonMembership)
G_NonMembership_Token=[ ];
for i=[1:Standard_Matrix_Row]
    for j=[1:Standard_Matrix_Column]
        G_NonMembership_Token(end+1)=(1/max(find(Reverse_Standard_Matrix{i,j}{1})))*sum((Reverse_Standard_Matrix{i,j}{1}).^0.1)*(1/0.1);
    end
end
G_NonMembership=[ ];
for i=[1:Standard_Matrix_Row:Standard_Matrix_Row*Standard_Matrix_Column]
    Row=[G_NonMembership_Token(i),G_NonMembership_Token(i+1),G_NonMembership_Token(i+2),G_NonMembership_Token(i+3)];
    G_NonMembership=[G_NonMembership;Row];
end
display(G_NonMembership)
S_NonMembership=V_NonMembership./(V_NonMembership+G_NonMembership)
E_NonMembership_Token=[ ];
for i=[1:Standard_Matrix_Row]
    for j=[1:Standard_Matrix_Column]
        if max(find(Reverse_Standard_Matrix{i,j}{1}))==1
            E_NonMembership_Token(end+1)=Reverse_Standard_Matrix{i,j}{1}(1);
        elseif  max(find(Reverse_Standard_Matrix{i,j}{1}))==2
            E_NonMembership_Token(end+1)=(1/2)*abs(Reverse_Standard_Matrix{i,j}{1}(2)-Reverse_Standard_Matrix{i,j}{1}(1));
        else if max(find(Reverse_Standard_Matrix{i,j}{1}))==3
                E_NonMembership_Token(end+1)=(1/3)*(abs(Reverse_Standard_Matrix{i,j}{1}(1)-Reverse_Standard_Matrix{i,j}{1}(2))+...
                    abs(Reverse_Standard_Matrix{i,j}{1}(2)-Reverse_Standard_Matrix{i,j}{1}(3)));
            end
        end
    end
end
E_NonMembership=[ ];
for i=[1:Standard_Matrix_Row:Standard_Matrix_Row*Standard_Matrix_Column]
    Row=[E_NonMembership_Token(i),E_NonMembership_Token(i+1),E_NonMembership_Token(i+2),E_NonMembership_Token(i+3)];
    E_NonMembership=[E_NonMembership;Row];
end
display(E_NonMembership)
S_Matrix=S_Membership+S_NonMembership
E_Matrix=E_Membership+E_NonMembership
Eta=0.5
% Eta is the preference of precept's score and the amount of its information.
Omega=[ ];
% Omega is the weight.
for j=[1:Standard_Matrix_Column]
    Omega(end+1)=(Eta*(sum(1-S_Matrix(:,j)))+(1-Eta)*(sum(1-E_Matrix(:,j))))/...
        ((Eta*(sum(1-S_Matrix(:,1)))+(1-Eta)*(sum(1-E_Matrix(:,1))))+(Eta*(sum(1-S_Matrix(:,2)))+(1-Eta)*(sum(1-E_Matrix(:,2))))+...
        (Eta*(sum(1-S_Matrix(:,3)))+(1-Eta)*(sum(1-E_Matrix(:,3))))+(Eta*(sum(1-S_Matrix(:,4)))+(1-Eta)*(sum(1-E_Matrix(:,4)))));
end
display(Omega)
Best_Membership_Token=[ ];
Worst_Membership_Token=[ ];
Best_NonMembership_Token=[ ];
Worst_NonMembership_Token=[ ];
for j=[1:Standard_Matrix_Column]
    for i=[1:Standard_Matrix_Row]
        Best_Membership_Token(end+1)=[max(Standard_Matrix{i,j}{1})];
        Worst_Membership_Token(end+1)=[min(Standard_Matrix{i,j}{1})];
        Best_NonMembership_Token(end+1)=[min(Standard_Matrix{i,j}{2})];
        Worst_NonMembership_Token(end+1)=[max(Standard_Matrix{i,j}{2})];
    end
end
Best_Membership=[ ];
Worst_Membership=[ ];
Best_NonMembership=[ ];
Worst_NonMembership=[ ];
Count_Number=1+Standard_Matrix_Column*(Standard_Matrix_Row-1);
for i=[1:Standard_Matrix_Column:Count_Number]
    row1=max(Best_Membership_Token(i:i+Standard_Matrix_Column-1));
    Best_Membership=[Best_Membership,row1];
    row2=min(Worst_Membership_Token(i:i+Standard_Matrix_Column-1));
    Worst_Membership=[Worst_Membership,row2];
    row3=min(Best_NonMembership_Token(i:i+Standard_Matrix_Column-1));
    Best_NonMembership=[Best_NonMembership,row3];
    row4=max(Worst_NonMembership_Token(i:i+Standard_Matrix_Column-1));
    Worst_NonMembership=[Worst_NonMembership,row4];
end
display(Best_Membership)
display(Worst_Membership)
display(Best_NonMembership)
display(Worst_NonMembership)
Distance_Positive_Part_1_Token=[ ];
Distance_Positive_Part_2_Token=[ ];
for j=[1:Standard_Matrix_Column]
    for i=[1:Standard_Matrix_Row]
        if max(find(Standard_Matrix{i,j}{1}))==1
            Distance_Positive_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Best_Membership(j))^2)/...
                max(find(Standard_Matrix{i,j}{1})))^0.5;
        else if max(find(Standard_Matrix{i,j}{1}))==2
                Distance_Positive_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Best_Membership(j))^2+...
                    (Standard_Matrix{i,j}{1}(2)-Best_Membership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{1})))^0.5;
            else if max(find(Standard_Matrix{i,j}{1}))==3
                    Distance_Positive_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Best_Membership(j))^2+...
                        (Standard_Matrix{i,j}{1}(2)-Best_Membership(j))^2+...
                        (Standard_Matrix{i,j}{1}(3)-Best_Membership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{1})))^0.5;
                end
            end
        end
    end
end
for j=[1:Standard_Matrix_Column]
    for i=[1:Standard_Matrix_Row]
        if max(find(Standard_Matrix{i,j}{2}))==1
            Distance_Positive_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Best_NonMembership(j))^2)/...
                max(find(Standard_Matrix{i,j}{2})))^0.5;
        else if max(find(Standard_Matrix{i,j}{2}))==2
                Distance_Positive_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Best_NonMembership(j))^2+...
                    (Standard_Matrix{i,j}{2}(2)-Best_NonMembership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{2})))^0.5;
            else if max(find(Standard_Matrix{i,j}{2}))==3
                    Distance_Positive_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Best_NonMembershi(j))^2+...
                        (Standard_Matrix{i,j}{2}(2)-Best_NonMembership(j))^2+...
                        (Standard_Matrix{i,j}{2}(3)-Best_NonMembership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{2})))^0.5;
                end
            end
        end
    end
end
Distance_Negative_Part_1_Token=[ ];
Distance_Negative_Part_2_Token=[ ];
for j=[1:Standard_Matrix_Column]
    for i=[1:Standard_Matrix_Row]
        if max(find(Standard_Matrix{i,j}{1}))==1
            Distance_Negative_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Worst_Membership(j))^2)/...
                max(find(Standard_Matrix{i,j}{1})))^0.5;
        else if max(find(Standard_Matrix{i,j}{1}))==2
                Distance_Negative_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Worst_Membership(j))^2+...
                    (Standard_Matrix{i,j}{1}(2)-Worst_Membership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{1})))^0.5;
            else if max(find(Standard_Matrix{i,j}{1}))==3
                    Distance_Negative_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Worst_Membership(j))^2+...
                        (Standard_Matrix{i,j}{1}(2)-Worst_Membership(j))^2+...
                        (Standard_Matrix{i,j}{1}(3)-Worst_Membership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{1})))^0.5;
                end
            end
        end
    end
end
for j=[1:Standard_Matrix_Column]
    for i=[1:Standard_Matrix_Row]
        if max(find(Standard_Matrix{i,j}{2}))==1
            Distance_Negative_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Worst_NonMembership(j))^2)/...
                max(find(Standard_Matrix{i,j}{2})))^0.5;
        else if max(find(Standard_Matrix{i,j}{2}))==2
                Distance_Negative_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Worst_NonMembership(j))^2+...
                    (Standard_Matrix{i,j}{2}(2)-Worst_NonMembership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{2})))^0.5;
            else if max(find(Standard_Matrix{i,j}{2}))==3
                    Distance_Negative_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Worst_NonMembership(j))^2+...
                        (Standard_Matrix{i,j}{2}(2)-Worst_NonMembership(j))^2+...
                        (Standard_Matrix{i,j}{2}(3)-Worst_NonMembership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{2})))^0.5;
                end
            end
        end
    end
end
Distance_Positive_Part_1=[ ];
Distance_Positive_Part_2=[ ];
for i=[1:Standard_Matrix_Column]
    row1=[Distance_Positive_Part_1_Token(i),...
        Distance_Positive_Part_1_Token(i+Standard_Matrix_Column),...
        Distance_Positive_Part_1_Token(i+2*Standard_Matrix_Column),...
        Distance_Positive_Part_1_Token(i+3*Standard_Matrix_Column)];
    row2=[Distance_Positive_Part_2_Token(i),...
        Distance_Positive_Part_2_Token(i+Standard_Matrix_Column),...
        Distance_Positive_Part_2_Token(i+2*Standard_Matrix_Column),...
        Distance_Positive_Part_2_Token(i+3*Standard_Matrix_Column)];
    Distance_Positive_Part_1=[Distance_Positive_Part_1;row1];
    Distance_Positive_Part_2=[Distance_Positive_Part_2;row1];
end
Distance_Positive=0.5*(Distance_Positive_Part_1+Distance_Positive_Part_2)
Distance_Negative_Part_1=[ ];
Distance_Negative_Part_2=[ ];
for i=[1:Standard_Matrix_Column]
    row1=[Distance_Negative_Part_1_Token(i),...
        Distance_Negative_Part_1_Token(i+Standard_Matrix_Column),...
        Distance_Negative_Part_1_Token(i+2*Standard_Matrix_Column),...
        Distance_Negative_Part_1_Token(i+3*Standard_Matrix_Column)];
    row2=[Distance_Negative_Part_2_Token(i),...
        Distance_Negative_Part_2_Token(i+Standard_Matrix_Column),...
        Distance_Negative_Part_2_Token(i+2*Standard_Matrix_Column),...
        Distance_Negative_Part_2_Token(i+3*Standard_Matrix_Column)];
    Distance_Negative_Part_1=[Distance_Negative_Part_1;row1];
    Distance_Negative_Part_2=[Distance_Negative_Part_2;row1];
end
Distance_Negative=0.5*(Distance_Negative_Part_1+Distance_Negative_Part_2)
Proximity_Matrix=Distance_Negative./(Distance_Positive+Distance_Negative)
Best_Proximity_Token=max(Proximity_Matrix);
Delta=0.3
% Delta is regret avoidance factor.
Best_Proximity=repmat(Best_Proximity_Token,Standard_Matrix_Row,1)
Regret_Perception_Matrix=Proximity_Matrix+1-exp(-Delta*(Proximity_Matrix-Best_Proximity))
Best_Regret_Perception=max(Regret_Perception_Matrix)
Worst_Regret_Perception=min(Regret_Perception_Matrix)
s_Token=[ ];
% s are parts of population utility value.
for j=[1:Standard_Matrix_Column]
    for i=[1:Standard_Matrix_Row]
        s_Token(end+1)=Omega(j)*(Best_Regret_Perception(j)-Regret_Perception_Matrix(i,j))./...
            (Best_Regret_Perception(j)-Worst_Regret_Perception(j));
    end
end
s=[ ];
for i=[1:Standard_Matrix_Row]
    row=[s_Token(i),s_Token(i+Standard_Matrix_Column),...
        s_Token(i+2*Standard_Matrix_Column),s_Token(i+3*Standard_Matrix_Column)];
    s=[s;row];
end
display(s)
S=[ ];
% S is population utility value.
for i=[1:Standard_Matrix_Row]
    S(end+1)=[sum(s(i,:))];
end
display(S)
R=[ ];
% R is individual regret value.
for i=[1:Standard_Matrix_Row]
    R(end+1)=[max(s(i,:))];
end
display(R)
Best_S=min(S)
Worst_S=max(S)
Best_R=min(R)
Worst_R=max(R)
v=0.5
% v is decision mechanism coefficient.
Q=[ ];
% Q is the evaluation index.
for i=[1:Standard_Matrix_Row]
    Q(end+1)=v*(S(i)-Best_S)/(Worst_S-Best_S)+(1-v)*(R(i)-Best_R)/(Worst_R-Best_R);
end
display(Q)
Rank_Q=sort(Q)
Q_4th=find(Q==max(Rank_Q));
Q_3rd=find(Q==Rank_Q(end-1));
Q_2nd=find(Q==Rank_Q(end-2));
Q_1st=find(Q==Rank_Q(end-3));
Rank=[Q_1st,Q_2nd,Q_3rd,Q_4th]
Result='The best plan is A';
Result=[Result,num2str(find(Q==min(Q)))];
display(Result)


% Following operations are sensitivity analysis 1. (v)


vQ=[ ];
for Disturbed_v_Token=[0:0.1:1]
    for i=[1:Standard_Matrix_Row]
        vQ(end+1)=Disturbed_v_Token*(S(i)-Best_S)/(Worst_S-Best_S)+...
            (1-Disturbed_v_Token)*(R(i)-Best_R)/(Worst_R-Best_R);
    end
end
Disturbed_v=[0:0.1:1];
Disturbed_v_Size=size(Disturbed_v);
Disturbed_v_Column=Disturbed_v_Size(2);
Sensitive_Analysis_Matrix_1=[ ];
Count_Number=1+Standard_Matrix_Column*(Disturbed_v_Column-1);
for i=[1:Standard_Matrix_Column:Count_Number]
    row=vQ(i:i+Standard_Matrix_Column-1);
    Sensitive_Analysis_Matrix_1=[Sensitive_Analysis_Matrix_1;row];
end
display(Sensitive_Analysis_Matrix_1)
Sensitive_Analysis_Matrix_1_Size=size(Sensitive_Analysis_Matrix_1);
Sensitive_Analysis_Matrix_1_Row=Sensitive_Analysis_Matrix_1_Size(1);
Sensitive_Analysis_Matrix_1_Column=Sensitive_Analysis_Matrix_1_Size(2);
Rank_Matrix_Q_1st=[ ];
Rank_Matrix_Q_2nd=[ ];
Rank_Matrix_Q_3rd=[ ];
Rank_Matrix_Q_4th=[ ];
for i=[1:Sensitive_Analysis_Matrix_1_Row]
    Each_Row=Sensitive_Analysis_Matrix_1(i,:);
    Rank_Row=sort(Sensitive_Analysis_Matrix_1(i,:));
    [row1,Q_4th]=find(Each_Row==max(Rank_Row));
    [row2,Q_3rd]=find(Each_Row==Rank_Row(end-1));
    [row3,Q_2nd]=find(Each_Row==Rank_Row(end-2));
    [row4,Q_1st]=find(Each_Row==Rank_Row(end-3));
    Rank_Matrix_Q_1st(end+1)=[Q_1st];
    Rank_Matrix_Q_2nd(end+1)=[Q_2nd];
    Rank_Matrix_Q_3rd(end+1)=[Q_3rd];
    Rank_Matrix_Q_4th(end+1)=[Q_4th];
end
Rank_Sensitive_Analysis_Matrix_1=[ ];
for i=[1:Sensitive_Analysis_Matrix_1_Row]
    row=[Rank_Matrix_Q_1st(i), Rank_Matrix_Q_2nd(i),Rank_Matrix_Q_3rd(i),Rank_Matrix_Q_4th(i)];
    Rank_Sensitive_Analysis_Matrix_1=[Rank_Sensitive_Analysis_Matrix_1;row];
end
display(Rank_Sensitive_Analysis_Matrix_1)


% Following operations are sensitivity analysis 2. (weight)


Xi=[0:0.05:2]
% Xi is disturbance parameters.
k_Token=[ ];
% k is the weight after disturbing.
Omega_Size=size(Omega);
Omega_Column=Omega_Size(2);
Xi_Size=size(Xi);
Xi_Row=Xi_Size(1);
Xi_Column=Xi_Size(2);
for i=[1:Omega_Column]
    for j=[1:Xi_Column]
        k_Token(end+1)=[(1-Omega(i)*Xi(j))/(1-Omega(i))];
    end
end
k=[ ];
Count_Number=1+Standard_Matrix_Column*(Standard_Matrix_Column-1);
for i=[1:Omega_Column:Count_Number]
    row=k_Token(i:i+Xi_Column-1);
    k=[k;row];
end
display(k)
k_Size=size(k);
k_Row=k_Size(1);
k_Column=k_Size(2);
Disturbed_Omega_Token=[ ];
n=1;
for i=[1:k_Row]
    for j=[1:k_Column]
        for m=[1:Omega_Column]
            if m==n
                Omega_Token=Omega(m)*Xi(j);
            else  Omega_Token=Omega(m)*k(i,j);
            end
            Disturbed_Omega_Token(end+1)=[Omega_Token];
        end
    end
    n=n+1;
end
Disturbed_Omega_Matrix=[ ];
Disturbed_Omega_Matrix_Size=size(Disturbed_Omega_Token);
Disturbed_Omega_Matrix_Column=Disturbed_Omega_Matrix_Size(2);
for i=[1:Omega_Column:Disturbed_Omega_Matrix_Column]
    row=Disturbed_Omega_Token(i:i+3);
    Disturbed_Omega_Matrix=[Disturbed_Omega_Matrix;row];
end
display(Disturbed_Omega_Matrix)
Original_Matrix=...
    {{0.2,[0.6,0.8]},{[0.6,0.8],[0.1,0.3]},{[0.4,0.5],[0.2,0.3]},{[0.5,0.6],0.5};...
    {[0.3,0.5],0.5},{[0.3,0.4,0.6],0.6},{0.2,[0.5,0.6]},{[0.5,0.7],[0.4,0.5]};...
    {0.7,[0.3,0.4,0.5]},{[0.5,0.6,0.7],0.7},{[0.7,0.8],[0.4,0.5]},{[0.3,0.4],0.4};...
    {[0.4,0.5],[0.2,0.4]},{0.6,[0.2,0.3]},{[0.3,0.5,0.6],0.5},{[0.7,0.8],[0.2,0.3]}};
Sensitivity_Analysis_Matrix_2=[ ];
for ii=[1:Disturbed_Omega_Matrix_Column/4]
    Disturbed_Omega=[Disturbed_Omega_Matrix(ii,:)];
    Standard_Matrix=Original_Matrix;
    Standard_Matrix_Size=size(Standard_Matrix);
    Standard_Matrix_Row=Standard_Matrix_Size(1);
    Standard_Matrix_Column=Standard_Matrix_Size(2);
    for i=[1:Standard_Matrix_Row]
        Exchange_Unit=Standard_Matrix{i,1}(1);
        Standard_Matrix{i,1}(1)=Standard_Matrix{i,1}(2);
        Standard_Matrix{i,1}(2)=Exchange_Unit;
    end
    Best_Membership_Token=[ ];
    Worst_Membership_Token=[ ];
    Best_NonMembership_Token=[ ];
    Worst_NonMembership_Token=[ ];
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            Best_Membership_Token(end+1)=[max(Standard_Matrix{i,j}{1})];
            Worst_Membership_Token(end+1)=[min(Standard_Matrix{i,j}{1})];
            Best_NonMembership_Token(end+1)=[min(Standard_Matrix{i,j}{2})];
            Worst_NonMembership_Token(end+1)=[max(Standard_Matrix{i,j}{2})];
        end
    end
    Best_Membership=[ ];
    Worst_Membership=[ ];
    Best_NonMembership=[ ];
    Worst_NonMembership=[ ];
    Count_Number=1+Standard_Matrix_Column*(Standard_Matrix_Row-1);
    for i=[1:4:Count_Number]
        row1=max(Best_Membership_Token(i:i+3));
        Best_Membership=[Best_Membership,row1];
        row2=min(Worst_Membership_Token(i:i+3));
        Worst_Membership=[Worst_Membership,row2];
        row3=min(Best_NonMembership_Token(i:i+3));
        Best_NonMembership=[Best_NonMembership,row3];
        row4=max(Worst_NonMembership_Token(i:i+3));
        Worst_NonMembership=[Worst_NonMembership,row4];
    end
    Distance_Positive_Part_1_Token=[ ];
    Distance_Positive_Part_2_Token=[ ];
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            if max(find(Standard_Matrix{i,j}{1}))==1
                Distance_Positive_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Best_Membership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{1})))^0.5;
            else if max(find(Standard_Matrix{i,j}{1}))==2
                    Distance_Positive_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Best_Membership(j))^2+...
                        (Standard_Matrix{i,j}{1}(2)-Best_Membership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{1})))^0.5;
                else if max(find(Standard_Matrix{i,j}{1}))==3
                        Distance_Positive_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Best_Membership(j))^2+...
                            (Standard_Matrix{i,j}{1}(2)-Best_Membership(j))^2+...
                            (Standard_Matrix{i,j}{1}(3)-Best_Membership(j))^2)/...
                            max(find(Standard_Matrix{i,j}{1})))^0.5;
                    end
                end
            end
        end
    end
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            if max(find(Standard_Matrix{i,j}{2}))==1
                Distance_Positive_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Best_NonMembership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{2})))^0.5;
            else if max(find(Standard_Matrix{i,j}{2}))==2
                    Distance_Positive_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Best_NonMembership(j))^2+...
                        (Standard_Matrix{i,j}{2}(2)-Best_NonMembership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{2})))^0.5;
                else if max(find(Standard_Matrix{i,j}{2}))==3
                        Distance_Positive_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Best_NonMembershi(j))^2+...
                            (Standard_Matrix{i,j}{2}(2)-Best_NonMembership(j))^2+...
                            (Standard_Matrix{i,j}{2}(3)-Best_NonMembership(j))^2)/...
                            max(find(Standard_Matrix{i,j}{2})))^0.5;
                    end
                end
            end
        end
    end
    Distance_Negative_Part_1_Token=[ ];
    Distance_Negative_Part_2_Token=[ ];
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            if max(find(Standard_Matrix{i,j}{1}))==1
                Distance_Negative_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Worst_Membership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{1})))^0.5;
            else if max(find(Standard_Matrix{i,j}{1}))==2
                    Distance_Negative_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Worst_Membership(j))^2+...
                        (Standard_Matrix{i,j}{1}(2)-Worst_Membership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{1})))^0.5;
                else if max(find(Standard_Matrix{i,j}{1}))==3
                        Distance_Negative_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Worst_Membership(j))^2+...
                            (Standard_Matrix{i,j}{1}(2)-Worst_Membership(j))^2+...
                            (Standard_Matrix{i,j}{1}(3)-Worst_Membership(j))^2)/...
                            max(find(Standard_Matrix{i,j}{1})))^0.5;
                    end
                end
            end
        end
    end
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            if max(find(Standard_Matrix{i,j}{2}))==1
                Distance_Negative_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Worst_NonMembership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{2})))^0.5;
            else if max(find(Standard_Matrix{i,j}{2}))==2
                    Distance_Negative_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Worst_NonMembership(j))^2+...
                        (Standard_Matrix{i,j}{2}(2)-Worst_NonMembership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{2})))^0.5;
                else if max(find(Standard_Matrix{i,j}{2}))==3
                        Distance_Negative_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Worst_NonMembership(j))^2+...
                            (Standard_Matrix{i,j}{2}(2)-Worst_NonMembership(j))^2+...
                            (Standard_Matrix{i,j}{2}(3)-Worst_NonMembership(j))^2)/...
                            max(find(Standard_Matrix{i,j}{2})))^0.5;
                    end
                end
            end
        end
    end
    Distance_Positive_Part_1=[ ];
    Distance_Positive_Part_2=[ ];
    for i=[1:Standard_Matrix_Column]
        row1=[Distance_Positive_Part_1_Token(i),...
            Distance_Positive_Part_1_Token(i+Standard_Matrix_Column),...
            Distance_Positive_Part_1_Token(i+2*Standard_Matrix_Column),...
            Distance_Positive_Part_1_Token(i+3*Standard_Matrix_Column)];
        row2=[Distance_Positive_Part_2_Token(i),...
            Distance_Positive_Part_2_Token(i+Standard_Matrix_Column),...
            Distance_Positive_Part_2_Token(i+2*Standard_Matrix_Column),...
            Distance_Positive_Part_2_Token(i+3*Standard_Matrix_Column)];
        Distance_Positive_Part_1=[Distance_Positive_Part_1;row1];
        Distance_Positive_Part_2=[Distance_Positive_Part_2;row1];
    end
    Distance_Positive=0.5*(Distance_Positive_Part_1+Distance_Positive_Part_2);
    Distance_Negative_Part_1=[ ];
    Distance_Negative_Part_2=[ ];
    for i=[1:Standard_Matrix_Column]
        row1=[Distance_Negative_Part_1_Token(i),...
            Distance_Negative_Part_1_Token(i+Standard_Matrix_Column),...
            Distance_Negative_Part_1_Token(i+2*Standard_Matrix_Column),...
            Distance_Negative_Part_1_Token(i+3*Standard_Matrix_Column)];
        row2=[Distance_Negative_Part_2_Token(i),...
            Distance_Negative_Part_2_Token(i+Standard_Matrix_Column),...
            Distance_Negative_Part_2_Token(i+2*Standard_Matrix_Column),...
            Distance_Negative_Part_2_Token(i+3*Standard_Matrix_Column)];
        Distance_Negative_Part_1=[Distance_Negative_Part_1;row1];
        Distance_Negative_Part_2=[Distance_Negative_Part_2;row1];
    end
    Distance_Negative=0.5*(Distance_Negative_Part_1+Distance_Negative_Part_2);
    Proximity_Matrix=Distance_Negative./(Distance_Positive+Distance_Negative);
    Best_Proximity=max(Proximity_Matrix);
    Delta=0.3;
    Regret_Perception_Matrix=Proximity_Matrix+1-exp(-Delta.*(Proximity_Matrix-Best_Proximity));
    Best_Regret_Perception=max(Regret_Perception_Matrix);
    Worst_Regret_Perception=min(Regret_Perception_Matrix);
    Regret_Perception_Matrix_Size=size(Regret_Perception_Matrix);
    Regret_Perception_Matrix_Row=Regret_Perception_Matrix_Size(1);
    Regret_Perception_Matrix_Column=Regret_Perception_Matrix_Size(2);
    s_Token=[ ];
    for j=[1:Regret_Perception_Matrix_Column]
        for i=[1:Regret_Perception_Matrix_Row]
            s_Token(end+1)=Disturbed_Omega(j)*(Best_Regret_Perception(j)-Regret_Perception_Matrix(i,j))./...
                (Best_Regret_Perception(j)-Worst_Regret_Perception(j));
        end
    end
    s=[ ];
    for i=[1:Regret_Perception_Matrix_Row]
        row=[s_Token(i),s_Token(i+4),s_Token(i+8),s_Token(i+12)];
        s=[s;row];
    end
    S=[ ];
    for i=[1:Regret_Perception_Matrix_Row]
        S(end+1)=[sum(s(i,:))];
    end
    R=[ ];
    for i=[1:Regret_Perception_Matrix_Row]
        R(end+1)=[max(s(i,:))];
    end
    Best_S=min(S);
    Worst_S=max(S);
    Best_R=min(R);
    Worst_R=max(R);
    v=0.5;
    Q=[ ];
    for i=[1:Regret_Perception_Matrix_Row]
        Q(end+1)=v*(S(i)-Best_S)/(Worst_S-Best_S)+(1-v)*(R(i)-Best_R)/(Worst_R-Best_R);
    end
    SAMrow=Q;
    Sensitivity_Analysis_Matrix_2=[Sensitivity_Analysis_Matrix_2;SAMrow];
end
display(Sensitivity_Analysis_Matrix_2)
Rank_Matrix_Q_1st=[ ];
Rank_Matrix_Q_2nd=[ ];
Rank_Matrix_Q_3rd=[ ];
Rank_Matrix_Q_4th=[ ];
Rank_Sensitivity_Analysis_Matrix_2=[ ];
for i=[1:ii]
    Each_Row=Sensitivity_Analysis_Matrix_2(i,:);
    Rank_Row=sort(Sensitivity_Analysis_Matrix_2(i,:));
    [row1,Q_4th]=find(Each_Row==max(Rank_Row));
    [row2,Q_3rd]=find(Each_Row==Rank_Row(end-1));
    [row3,Q_2nd]=find(Each_Row==Rank_Row(end-2));
    [row4,Q_1st]=find(Each_Row==Rank_Row(end-3));
    Rank_Matrix_Q_1st(end+1)=[Q_1st];
    Rank_Matrix_Q_2nd(end+1)=[Q_2nd];
    Rank_Matrix_Q_3rd(end+1)=[Q_3rd];
    Rank_Matrix_Q_4th(end+1)=[Q_4th];
end
for i=[1:ii]
    row=[Rank_Matrix_Q_1st(i), Rank_Matrix_Q_2nd(i),Rank_Matrix_Q_3rd(i),Rank_Matrix_Q_4th(i)];
    Rank_Sensitivity_Analysis_Matrix_2=[Rank_Sensitivity_Analysis_Matrix_2;row];
end
display(Rank_Sensitivity_Analysis_Matrix_2)
% The following steps aim to find the amount of Q1 to Q4 in column 1.
Q_1st_Amount_1=length(find(Rank_Sensitivity_Analysis_Matrix_2(:,1)==Rank_Sensitive_Analysis_Matrix_1((Disturbed_v_Column-1)/2,1)))
Q_2nd_Amount_1=length(find(Rank_Sensitivity_Analysis_Matrix_2(:,1)==Rank_Sensitive_Analysis_Matrix_1((Disturbed_v_Column-1)/2,2)))
Q_3rd_Amount_1=length(find(Rank_Sensitivity_Analysis_Matrix_2(:,1)==Rank_Sensitive_Analysis_Matrix_1((Disturbed_v_Column-1)/2,3)))
Q_4th_Amount_1=length(find(Rank_Sensitivity_Analysis_Matrix_2(:,1)==Rank_Sensitive_Analysis_Matrix_1((Disturbed_v_Column-1)/2,4)))
Q_1st_Amount_2=length(find(Rank_Sensitivity_Analysis_Matrix_2(:,2)==Rank_Sensitive_Analysis_Matrix_1((Disturbed_v_Column-1)/2,1)))
Q_2nd_Amount_2=length(find(Rank_Sensitivity_Analysis_Matrix_2(:,2)==Rank_Sensitive_Analysis_Matrix_1((Disturbed_v_Column-1)/2,2)))


% Following operations are sensitivity analysis 3. (regret avoidance factor)


R_Q=[ ];
R_Rank=[ ];
for Delta=[0:0.05:0.4]
    Regret_Perception_Matrix=Proximity_Matrix+1-exp(-Delta*(Proximity_Matrix-Best_Proximity));
    Best_Regret_Perception=max(Regret_Perception_Matrix);
    Worst_Regret_Perception=min(Regret_Perception_Matrix);
    s_Token=[ ];
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            s_Token(end+1)=Omega(j)*(Best_Regret_Perception(j)-Regret_Perception_Matrix(i,j))./...
                (Best_Regret_Perception(j)-Worst_Regret_Perception(j));
        end
    end
    s=[ ];
    for i=[1:Standard_Matrix_Row]
        row=[s_Token(i),s_Token(i+Standard_Matrix_Column),...
            s_Token(i+2*Standard_Matrix_Column),s_Token(i+3*Standard_Matrix_Column)];
        s=[s;row];
    end
    S=[ ];
    for i=[1:Standard_Matrix_Row]
        S(end+1)=[sum(s(i,:))];
    end
    R=[ ];
    for i=[1:Standard_Matrix_Row]
        R(end+1)=[max(s(i,:))];
    end
    Best_S=min(S);
    Worst_S=max(S);
    Best_R=min(R);
    Worst_R=max(R);
    v=0.5;
    Q=[ ];
    for i=[1:Standard_Matrix_Row]
        Q(end+1)=v*(S(i)-Best_S)/(Worst_S-Best_S)+(1-v)*(R(i)-Best_R)/(Worst_R-Best_R);
    end
    Rank_Q=sort(Q);
    Q_4th=find(Q==max(Rank_Q));
    Q_3rd=find(Q==Rank_Q(end-1));
    Q_2nd=find(Q==Rank_Q(end-2));
    Q_1st=find(Q==Rank_Q(end-3));
    Rank=[Q_1st,Q_2nd,Q_3rd,Q_4th];
    Q_Change=Q;
    Rank_Change=Rank;
    R_Q=[R_Q;Q_Change];
    R_Rank=[R_Rank;Rank_Change];
end
display(R_Q)
display(R_Rank)


% Following operations are sensitivity analysis 4. (weight avoidance factor)


Omega=[ ];
for Eta=[0:0.1:1]
    for j=[1:Standard_Matrix_Column]
        Omega(end+1)=(Eta*(sum(1-S_Matrix(:,j)))+(1-Eta)*(sum(1-E_Matrix(:,j))))/...
            ((Eta*(sum(1-S_Matrix(:,1)))+(1-Eta)*(sum(1-E_Matrix(:,1))))+(Eta*(sum(1-S_Matrix(:,2)))+(1-Eta)*(sum(1-E_Matrix(:,2))))+...
            (Eta*(sum(1-S_Matrix(:,3)))+(1-Eta)*(sum(1-E_Matrix(:,3))))+(Eta*(sum(1-S_Matrix(:,4)))+(1-Eta)*(sum(1-E_Matrix(:,4)))));
    end
end
Eta=[0:0.1:1];
Size_Eta=size(Eta);
Eta_Column=Size_Eta(2);
Omega_Change=[ ];
Omega_Avoidance=[ ];
for i=[1:Omega_Column:Omega_Column*Eta_Column]
    Omega_Change=Omega(i:i+Omega_Column-1);
    Omega_Avoidance=[Omega_Avoidance;Omega_Change];
end
display(Omega_Avoidance)
Omega=[ ];
Sensitivity_Analysis_Matrix_4=[ ];
Q_4=[ ];
Rank_4=[ ];
for ii=[1:Eta_Column]
    Omega=Omega_Avoidance(ii,:);
    Best_Membership_Token=[ ];
    Worst_Membership_Token=[ ];
    Best_NonMembership_Token=[ ];
    Worst_NonMembership_Token=[ ];
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            Best_Membership_Token(end+1)=[max(Standard_Matrix{i,j}{1})];
            Worst_Membership_Token(end+1)=[min(Standard_Matrix{i,j}{1})];
            Best_NonMembership_Token(end+1)=[min(Standard_Matrix{i,j}{2})];
            Worst_NonMembership_Token(end+1)=[max(Standard_Matrix{i,j}{2})];
        end
    end
    Best_Membership=[ ];
    Worst_Membership=[ ];
    Best_NonMembership=[ ];
    Worst_NonMembership=[ ];
    Count_Number=1+Standard_Matrix_Column*(Standard_Matrix_Row-1);
    for i=[1:Standard_Matrix_Column:Count_Number]
        row1=max(Best_Membership_Token(i:i+Standard_Matrix_Column-1));
        Best_Membership=[Best_Membership,row1];
        row2=min(Worst_Membership_Token(i:i+Standard_Matrix_Column-1));
        Worst_Membership=[Worst_Membership,row2];
        row3=min(Best_NonMembership_Token(i:i+Standard_Matrix_Column-1));
        Best_NonMembership=[Best_NonMembership,row3];
        row4=max(Worst_NonMembership_Token(i:i+Standard_Matrix_Column-1));
        Worst_NonMembership=[Worst_NonMembership,row4];
    end
    Distance_Positive_Part_1_Token=[ ];
    Distance_Positive_Part_2_Token=[ ];
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            if max(find(Standard_Matrix{i,j}{1}))==1
                Distance_Positive_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Best_Membership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{1})))^0.5;
            else if max(find(Standard_Matrix{i,j}{1}))==2
                    Distance_Positive_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Best_Membership(j))^2+...
                        (Standard_Matrix{i,j}{1}(2)-Best_Membership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{1})))^0.5;
                else if max(find(Standard_Matrix{i,j}{1}))==3
                        Distance_Positive_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Best_Membership(j))^2+...
                            (Standard_Matrix{i,j}{1}(2)-Best_Membership(j))^2+...
                            (Standard_Matrix{i,j}{1}(3)-Best_Membership(j))^2)/...
                            max(find(Standard_Matrix{i,j}{1})))^0.5;
                    end
                end
            end
        end
    end
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            if max(find(Standard_Matrix{i,j}{2}))==1
                Distance_Positive_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Best_NonMembership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{2})))^0.5;
            else if max(find(Standard_Matrix{i,j}{2}))==2
                    Distance_Positive_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Best_NonMembership(j))^2+...
                        (Standard_Matrix{i,j}{2}(2)-Best_NonMembership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{2})))^0.5;
                else if max(find(Standard_Matrix{i,j}{2}))==3
                        Distance_Positive_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Best_NonMembershi(j))^2+...
                            (Standard_Matrix{i,j}{2}(2)-Best_NonMembership(j))^2+...
                            (Standard_Matrix{i,j}{2}(3)-Best_NonMembership(j))^2)/...
                            max(find(Standard_Matrix{i,j}{2})))^0.5;
                    end
                end
            end
        end
    end
    Distance_Negative_Part_1_Token=[ ];
    Distance_Negative_Part_2_Token=[ ];
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            if max(find(Standard_Matrix{i,j}{1}))==1
                Distance_Negative_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Worst_Membership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{1})))^0.5;
            else if max(find(Standard_Matrix{i,j}{1}))==2
                    Distance_Negative_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Worst_Membership(j))^2+...
                        (Standard_Matrix{i,j}{1}(2)-Worst_Membership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{1})))^0.5;
                else if max(find(Standard_Matrix{i,j}{1}))==3
                        Distance_Negative_Part_1_Token(end+1)=(((Standard_Matrix{i,j}{1}(1)-Worst_Membership(j))^2+...
                            (Standard_Matrix{i,j}{1}(2)-Worst_Membership(j))^2+...
                            (Standard_Matrix{i,j}{1}(3)-Worst_Membership(j))^2)/...
                            max(find(Standard_Matrix{i,j}{1})))^0.5;
                    end
                end
            end
        end
    end
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            if max(find(Standard_Matrix{i,j}{2}))==1
                Distance_Negative_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Worst_NonMembership(j))^2)/...
                    max(find(Standard_Matrix{i,j}{2})))^0.5;
            else if max(find(Standard_Matrix{i,j}{2}))==2
                    Distance_Negative_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Worst_NonMembership(j))^2+...
                        (Standard_Matrix{i,j}{2}(2)-Worst_NonMembership(j))^2)/...
                        max(find(Standard_Matrix{i,j}{2})))^0.5;
                else if max(find(Standard_Matrix{i,j}{2}))==3
                        Distance_Negative_Part_2_Token(end+1)=(((Standard_Matrix{i,j}{2}(1)-Worst_NonMembership(j))^2+...
                            (Standard_Matrix{i,j}{2}(2)-Worst_NonMembership(j))^2+...
                            (Standard_Matrix{i,j}{2}(3)-Worst_NonMembership(j))^2)/...
                            max(find(Standard_Matrix{i,j}{2})))^0.5;
                    end
                end
            end
        end
    end
    Distance_Positive_Part_1=[ ];
    Distance_Positive_Part_2=[ ];
    for i=[1:Standard_Matrix_Column]
        row1=[Distance_Positive_Part_1_Token(i),...
            Distance_Positive_Part_1_Token(i+Standard_Matrix_Column),...
            Distance_Positive_Part_1_Token(i+2*Standard_Matrix_Column),...
            Distance_Positive_Part_1_Token(i+3*Standard_Matrix_Column)];
        row2=[Distance_Positive_Part_2_Token(i),...
            Distance_Positive_Part_2_Token(i+Standard_Matrix_Column),...
            Distance_Positive_Part_2_Token(i+2*Standard_Matrix_Column),...
            Distance_Positive_Part_2_Token(i+3*Standard_Matrix_Column)];
        Distance_Positive_Part_1=[Distance_Positive_Part_1;row1];
        Distance_Positive_Part_2=[Distance_Positive_Part_2;row1];
    end
    Distance_Positive=0.5*(Distance_Positive_Part_1+Distance_Positive_Part_2);
    Distance_Negative_Part_1=[ ];
    Distance_Negative_Part_2=[ ];
    for i=[1:Standard_Matrix_Column]
        row1=[Distance_Negative_Part_1_Token(i),...
            Distance_Negative_Part_1_Token(i+Standard_Matrix_Column),...
            Distance_Negative_Part_1_Token(i+2*Standard_Matrix_Column),...
            Distance_Negative_Part_1_Token(i+3*Standard_Matrix_Column)];
        row2=[Distance_Negative_Part_2_Token(i),...
            Distance_Negative_Part_2_Token(i+Standard_Matrix_Column),...
            Distance_Negative_Part_2_Token(i+2*Standard_Matrix_Column),...
            Distance_Negative_Part_2_Token(i+3*Standard_Matrix_Column)];
        Distance_Negative_Part_1=[Distance_Negative_Part_1;row1];
        Distance_Negative_Part_2=[Distance_Negative_Part_2;row1];
    end
    Distance_Negative=0.5*(Distance_Negative_Part_1+Distance_Negative_Part_2);
    Proximity_Matrix=Distance_Negative./(Distance_Positive+Distance_Negative);
    Best_Proximity=max(Proximity_Matrix);
    Delta=0.3;
    Regret_Perception_Matrix=Proximity_Matrix+1-exp(-Delta*(Proximity_Matrix-Best_Proximity));
    Best_Regret_Perception=max(Regret_Perception_Matrix);
    Worst_Regret_Perception=min(Regret_Perception_Matrix);
    s_Token=[ ];
    for j=[1:Standard_Matrix_Column]
        for i=[1:Standard_Matrix_Row]
            s_Token(end+1)=Omega(j)*(Best_Regret_Perception(j)-Regret_Perception_Matrix(i,j))./...
                (Best_Regret_Perception(j)-Worst_Regret_Perception(j));
        end
    end
    s=[ ];
    for i=[1:Standard_Matrix_Row]
        row=[s_Token(i),s_Token(i+Standard_Matrix_Column),...
            s_Token(i+2*Standard_Matrix_Column),s_Token(i+3*Standard_Matrix_Column)];
        s=[s;row];
    end
    S=[ ];
    for i=[1:Standard_Matrix_Row]
        S(end+1)=[sum(s(i,:))];
    end
    R=[ ];
    for i=[1:Standard_Matrix_Row]
        R(end+1)=[max(s(i,:))];
    end
    Best_S=min(S);
    Worst_S=max(S);
    Best_R=min(R);
    Worst_R=max(R);
    v=0.5;
    Q=[ ];
    for i=[1:Standard_Matrix_Row]
        Q(end+1)=v*(S(i)-Best_S)/(Worst_S-Best_S)+(1-v)*(R(i)-Best_R)/(Worst_R-Best_R);
    end
    Rank_Q=sort(Q);
    Q_4th=find(Q==max(Rank_Q));
    Q_3rd=find(Q==Rank_Q(end-1));
    Q_2nd=find(Q==Rank_Q(end-2));
    Q_1st=find(Q==Rank_Q(end-3));
    Rank=[Q_1st,Q_2nd,Q_3rd,Q_4th];
    Q_4=[Q_4;Q];
    Rank_4=[Rank_4;Rank];
end
display(Q_4)
display(Rank_4)