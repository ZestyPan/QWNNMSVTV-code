% function [ EPat W ] = PatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par )
% 
%         EPat = zeros(size(CurPat));
%         W    = zeros(size(CurPat));
%         for  i      =  1 : length(Self_arr)                                 % For each keypatch group
%             Temp    =   CurPat(:, NL_mat(1:Par.patnum,i));                  % Non-local similar patches to the keypatch
%             M_Temp  =   repmat(mean( Temp, 2 ),1,Par.patnum);
%             Temp    =   Temp-M_Temp;
% 
%             E_Temp 	=   WNNM( Temp, Par.c, Sigma_arr(Self_arr(i)), M_Temp, Par.ReWeiIter); % WNNM Estimation
%             EPat(:,NL_mat(1:Par.patnum,i))  = EPat(:,NL_mat(1:Par.patnum,i))+E_Temp;      
%             W(:,NL_mat(1:Par.patnum,i))     = W(:,NL_mat(1:Par.patnum,i))+ones(Par.patsize*Par.patsize,size(NL_mat(1:Par.patnum,i),1));
%         end
% end

function [EPat, W] = PatEstimation(NL_mat, Self_arr, Sigma_arr, CurPat, Par,E_Img)

    % 初始化估计补丁和权重矩阵
    EPat = zeros(size(CurPat));
    W    = zeros(size(CurPat));
    
    % 循环处理每个关键补丁组
    for i = 1:length(Self_arr)
        % 获取非局部相似补丁
        Temp = CurPat(:, NL_mat(1:Par.patnum, i));
        
        % 计算均值
        M_Temp = repmat(mean(Temp, 2), 1, Par.patnum);
        Temp = Temp - M_Temp;
        
        % WNNM 估计
        E_Temp = WNNM(Temp, Par.c, Sigma_arr(Self_arr(i)), M_Temp, Par.ReWeiIter,E_Img);
        
        % 更新 EPat
        EPat(:, NL_mat(1:Par.patnum, i)) = EPat(:, NL_mat(1:Par.patnum, i)) + E_Temp;
        
        % 更新 W
        % 确保 ones 的尺寸与要更新的 W 的部分相同
        W(:, NL_mat(1:Par.patnum, i)) = W(:, NL_mat(1:Par.patnum, i)) + ones(size(E_Temp, 1), 1);
    end
end