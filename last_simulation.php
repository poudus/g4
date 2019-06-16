<?php
	$connstring = 'dbname=g4 user=g4';
	$conn = pg_connect($connstring) or die ("could not connect");

	$query_max = 'select max(id) from g4s.simulation';
	$max_rs = pg_query($conn, $query_max);
	while ($rw = pg_fetch_row($max_rs)) {
		$sid = $rw[0];
	}

	$query = 'select * from g4s.frame where sid = ' . $sid . ' order by epoch, pid asc';
	$rs = pg_query($conn, $query);

$frames = array();

	while ($row = pg_fetch_row($rs)) {
		$frames[] = $row;
	}
	echo json_encode($frames);
	pg_close($conn);
	http_response_code(200);
?>
