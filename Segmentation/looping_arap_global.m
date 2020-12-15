% Todas as funcoes (exceto print3Dopt_grid) utilizadas aqui s?o do GPTOOLBOX 
% ARAP por segmento sem considerar pontos da superficie que precisam de
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


for h = 1:10
    
    % % % % Criterio de parada
    if norm(Xmin) < 1e-4
        break;
    end
    
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
        [~,Vnew,~,~] = print3Dopt_grid(V,F(seg==k,:),'zmin',...
            min(V(:,3)),'theta1_max',pi/180,'theta2_max',pi/180);
        toc;
        % % % % %
        
        segs{k} = Vnew; % Armazena os pontos da superficie apos metodo glo-
                        % bal aplicado aa segmenta??o

        b=unique([F(seg==k,1);F(seg==k,2);F(seg==k,3)],'rows'); % Tomando os indices do segmento

        bc=segs{k}(b,:); % Tomando os vertices do segmento apos rotacao global do segmento.
            
        % Aplicando ARAP ao segmento
        % % % %
        % % % % A partir dos vetores overhang_b e overhang_bc, em que
        % % % % indicam os indices dos pontos e as posicoes finais dos
        % % % % pontos que queremos transladar, respectivamente.
        % % % %
        fprintf('ARAP: ');
        tic;
        U = arap(V,F,b,bc,'Energy','elements'); % Processo Arap
        toc;
        V=U; % Atualiza??o
        fprintf('---------------------------');
    end
end
fprintf('\n');
tsurf(F,V);axis equal;view([1 0 0]);title('Global+ARAP'); %Imprimir superficie apos arap
figure;
tsurf(F,V_global);axis equal;view([1 0 0]);title('Global');%Imprimir superficie apos rotacao global.