% Todas as funcoes (exceto print3Dopt_grid e normalsurf) utilizadas aqui s?o do GPTOOLBOX 
% ARAP por segmento considerando apenas pontos da superficie que precisam de
% suporte.
clear;clc;
[V,F] = load_mesh('../MeshsegBenchmark-1.0/data/off/200.off'); % Carrega superficie em V e F
V_original=V; % Armazena os pontos da superficie na posicao original
[Xmin,V,~,~] = print3Dopt_grid(V,F,'zmin',...
        min(V(:,3))); % Aplica o metodo de rotacao global

% Escolhe a segmentacao da superficie e a imprime.
seg = plot_mesh_segmentation('../MeshsegBenchmark-1.0/data/off/200.off','../MeshsegBenchmark-1.0/data/seg/Benchmark/200/200_0.seg');
hold off;
% % % % % 
V_global=V; % Armazena os pontos da superf?cie apos sofrer rotacao global

bool_stop=0; % Variavel do crit?rio de parada

for h = 1:10
    
    if h ~= 1
        [Xmin,V,~,~] = print3Dopt_grid(V,F,'zmin',...
            min(V(:,3))); % Aplica o metodo da rota??o global
    end
    
    fprintf('\n__________________________________________________\n\n')
    fprintf(' Fase %d: Angulos de variacao: %.4f, %.4f',h,Xmin(1),Xmin(2));
    fprintf('\n__________________________________________________\n')
    
% % % % Inicio do processo de rota??o por segmentos
    for k = 1:max(seg)
        fprintf('\nFase %d: Carregando %d of %d. \n',h, k,max(seg))
        fprintf('Global: ');
        
        tic; 
        % Aplica o metodo de rotacao global no segmento
        [Xmin,Vnew,~,~] = print3Dopt_grid(V,F(seg==k,:),'zmin',...
            min(V(:,3)),'theta1_max',pi/90,'theta2_max',pi/90);
        % % % % %
        toc;
        
        segs{k} = Vnew; % Armazena os pontos da superficie apos metodo glo-
                        % bal aplicado aa segmenta??o
                        
        % Avaliacao para criterio de parada:
        % % % %
        % % % % Se os angulos de rotacao de cada segmento forem pequenos o
        % % % % suficiente, fazemos bool_stop = 1. Caso contrario, isto eh,
        % % % % se um dos dois angulos nao for pequeno o suficiente, a
        % % % % variavel bool_stop mantem-se igual a zero.
        % % % %
        if abs(Xmin(1)) < 1e-2 && abs(Xmin(2)) < 1e-2
            bool_stop = 1;
        else
            bool_stop = 0;
        end
        % % % % Fim Avaliacao para criterio de parada % % % %
        
        % % % % Tomando os indices do segmento
        b=unique([F(seg==k,1);F(seg==k,2);F(seg==k,3)],'rows');
        % % % %
        
        % Tomando os vertices do segmento apos rotacao global do segmento.
        bc=segs{k}(b,:);
        % % % %
        
        % % % % Avaliando a superficie quanto aa necessidade de suporte% %%
                                                                          %
        N_bc = normalsurf(Vnew,F); % Encontrando campo normal aa superficie
        N_bc = N_bc(b,:); % Tomando o campo normal apenas do seguemento.  %
                                                                          %
        % % % % Tomando a proje??o do campo normal do segmento.           %
        proj_N_bc = [N_bc(:,1) N_bc(:,2) zeros(length(N_bc),1)];          %
        % % % %                                                           %
                                                                          %
        alpha=acos(normrow(proj_N_bc)); % Calculando o angulo da normal   %
                                                                          %
        % Verificando se a normal esta apontada para cima ou para baixo:  %
        % % % %                                                           %
        % % % % Se a normal estiver apontada para baixo, poderemos        %
        % % % % ter a formacao de suportes externos. Caso a normal esteja %
        % % % % apontada para cima, poderemos ter a forma??o de suportes  %
        % % % % internos.                                                 %
        % % % %                                                           %
        v=-(N_bc(:,3) < 0);                                               %
        u=N_bc(:,3) > 0;                                                  %
        v=v+u;                                                            %
        alpha=v.*alpha; % Mudando o sinal do angulo da normal.            %
        % % % %                                                           %
                                                                          %
        bolean_overhang=alpha < -pi/4; % Verificando a necessidade de     %
                                       % suporte                          %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %               
        
        % % % %
        % % % % Criando o vetor de indices (overhang_b) dos pontos que
        % % % % queremos transladar e o vetor de posicao final dos pontos 
        % % % % (overhang_bc) que serao transladados pelo ARAP.
        % % % %
        overhang_b = zeros(sum(bolean_overhang),1);
        overhang_bc = zeros(sum(bolean_overhang),3);
        j=1;
        for i = 1:length(b)
            if bolean_overhang(i) == 1
                overhang_b(j) = b(i);
                overhang_bc(j,:) = segs{k}(b(i),:);
                j=j+1;
            end
        end
        % % % %
        
        % Aplicando ARAP ao segmento
        % % % %
        % % % % A partir dos vetores overhang_b e overhang_bc, em que
        % % % % indicam os indices dos pontos e as posicoes finais dos
        % % % % pontos que queremos transladar, respectivamente.
        % % % %
        fprintf('ARAP: ');
        tic;
        U = arap(V,F,overhang_b,overhang_bc,'Energy','elements'); % Processo Arap
        toc;
        V=U; % Atualizacao dos pontos
        fprintf('---------------------------');
    end
    
    % Criterio de parada
    % % % %
    % % % % Se ao final de cada analise dos segmentos, tiver pelo menos uma
    % % % % mudanca na rotacao global dos segmentos, a variavel bool_stop
    % % % % tera valor reinicializado (ou seja, 0). Se ao passarmos por
    % % % % todos os segmentos a variavel bool_stop tiver valor 1, entao o
    % % % % programa para. Isso faz com que o programa n?o fique rodando
    % % % % ate o final sem necessidade, ja que nao vai estar fazendo
    % % % % nenhuma operacao de arap nem de rotacao global.
    % % % %
    if bool_stop == 1
        break;
    end
    
end
fprintf('\n');
tsurf(F,V);axis equal;view([1 0 0]);title('Original');
figure;
tsurf(F,V_global);axis equal;view([1 0 0]);title('Global');%Imprimir superficie apos rotacao global.