function u = SetupDoms( DomainCells,dom_op,xg,xc )

    Ng = length(xg);

    Nc = length(xc);

    NumInterfs = DomainCells{dom_op+1}.m_NumInfefaces;

    u= cell(2,1);

    DGauss = zeros(6*NumInterfs,Ng);

    DCheb = zeros(6*NumInterfs,Nc);

    XG =zeros(NumInterfs,Ng);

    YG =zeros(NumInterfs,Ng);

    NxG =zeros(NumInterfs,Ng);

    NyG =zeros(NumInterfs,Ng);

    JG =zeros(NumInterfs,Ng);

    JpG =zeros(NumInterfs,Ng);

    XC =zeros(NumInterfs,Nc);

    YC =zeros(NumInterfs,Nc);

    NxC =zeros(NumInterfs,Nc);

    NyC =zeros(NumInterfs,Nc);

    JC=zeros(NumInterfs,Nc);

    JpC =zeros(NumInterfs,Nc);

    for i=1:NumInterfs

        [x,y] = DomainCells{dom_op+1}.geo(xg,i);

        XG(i,:)= x;

        YG(i,:)= y;

        [x,y] = DomainCells{dom_op+1}.normal(xg,i);

        NxG(i,:)= x;

        NyG(i,:)= y;

        JG(i,:) = DomainCells{dom_op+1}.J(xg,i);

        JpG(i,:) = DomainCells{dom_op+1}.Jp(xg,i);

        [x,y] = DomainCells{dom_op+1}.geo(xc,i);

        XC(i,:)= x;

        YC(i,:)= y;

        [x,y] = DomainCells{dom_op+1}.normal(xc,i);

        NxC(i,:)= x;

        NyC(i,:)= y;

        JC(i,:) = DomainCells{dom_op+1}.J(xc,i);

        JpC(i,:) = DomainCells{dom_op+1}.Jp(xc,i);

    end

    DGauss=[XG;YG;NxG;NyG;JG;JpG];

    DCheb =[XC;YC;NxC;NyC;JC;JpC];

    u{1} = DGauss;

    u{2} = DCheb;



end
