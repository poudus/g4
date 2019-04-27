<?php
$sid = $_GET["sid"];

	$connstring = 'dbname=g4 user=g4';
	$conn = pg_connect($connstring) or die ("could not connect");
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
