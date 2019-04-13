
create schema g4;

create sequence g4.simulation_sequo start 1;

create table g4.simulation (
id int primary key,
ts bigint not null,
epochs int not null,
frames int not null,
gravity numeric(6,2) not null,
rebound numeric(6,2) not null,
transfer numeric(6,2) not null,
views int not null,
likes int not null
);

create table g4.frame (
sid int not null,
epoch int not null,
x int not null,
y int not null,
radius numeric(6,2) not null,
color char(8) not null
);

