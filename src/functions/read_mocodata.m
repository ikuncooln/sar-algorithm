function [ MOCO_UNIT_HEAD, MOCO_UNIT ] = read_mocodata( file_name, pulse_count )
%READ_MOCODATA 读取运动补偿数据
%   补偿数据文件的存储格式为：MOCO_UNIT_HEAD（仅一个）, MOCO_UNIT（每个脉冲对应一个）
%
%   file_name mocodata文件名
%   pulse_count 要读取信息的脉冲个数，介于1到总的脉冲个数之前
%   MOCO_UNIT_HEAD struct 对应运补数据文件头
%   MOCO_UNIT struct 对应运补数据每个脉冲发射时位置信息
%   例：
%       [ MOCO_UNIT_HEAD, MOCO_UNIT ] = read_mocodata( file_name, 20480 );
%       disp(['载机参考prf', num2str(MOCO_UNIT_HEAD.ref_prf)]);
%       % 载机飞行空间轨迹
%       figure;
%       plot3(MOCO_UNIT.forward, MOCO_UNIT.cross, MOCO_UNIT.height);
%       xlabel('forward'); ylabel('cross'); zlabel('height');
%       grid on;

% 读取运动补偿数据
len_moco_unit_head = 9;
len_moco_unit = 12;
fid = fopen(file_name, 'r');
if fid==-1
    error(['文件', file_name, '不存在']);
end

d1 = fread(fid, len_moco_unit_head,'double');
MOCO_UNIT_HEAD = struct(...
    'origin_lon', d1(1),...
    'origin_lat', d1(2),...
    'scene_height', d1(3),...
    'ref_vel', d1(4),...
    'ref_prf', d1(5),...
    'baselength', d1(6),...
    'baseazlength', d1(7),...
    'baseangle', d1(8),...
    'anglefromeast', d1(9)...
);

d2 = fread(fid, [len_moco_unit, pulse_count],'double');
MOCO_UNIT = struct(...
    'time', d2(1, :),...
    'ref_lon', d2(2, :),...
    'ref_lat', d2(3, :),...
    'forward', d2(4, :),...
    'cross', d2(5, :),...
    'height', d2(6, :),...
    'ref_forward', d2(7, :),...
    'ref_east', d2(8, :),...
    'ref_north', d2(9, :),...
    'ref_height', d2(10, :),...
    'yaw', d2(11, :),...
    'pitch', d2(12, :)...
);

fclose(fid);

end

