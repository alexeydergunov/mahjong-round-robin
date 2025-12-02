# https://demonstrations.wolfram.com/SocialGolferProblem/
# https://www.mathpuzzle.com/MAA/54-Golf%20Tournaments/mathgames_08_14_07.html

def main():
    lines: list[str] = []
    with open("social_golfer_problem.txt", "r") as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                lines.append(line)

    player_index_map: dict[str, int] = {}
    line_index = 0
    while line_index < len(lines):
        line = lines[line_index]
        line_index += 1
        if line.endswith(":"):
            arr = line.split()
            player_count = int(arr[0])
            round_count = int(arr[2])
            player_index_map.clear()
            table_count = player_count // 4
            print(f"{player_count} players, {round_count} rounds:")
            player_opponents_map: dict[int, set[int]] = {}
            for player_index in range(player_count):
                player_opponents_map[player_index] = set()
            for round_index in range(round_count):
                line = lines[line_index + round_index]
                arr = line.split()
                assert len(arr) == table_count + 1
                assert arr[0] == "R" + str(round_index + 1) + ":"
                print(f"Round {round_index + 1}:")
                for table_index in range(1, table_count + 1):
                    print(f"  Table {table_index}:", end=" ")
                    for i in range(4):
                        c = arr[table_index][i]
                        if c not in player_index_map:
                            player_index = len(player_index_map)
                            player_index_map[c] = player_index
                        print(player_index_map[c] + 1, end=" ")
                    for i in range(4):
                        for j in range(4):
                            if i == j:
                                continue
                            player_index = player_index_map[arr[table_index][i]]
                            opponent_index = player_index_map[arr[table_index][j]]
                            player_opponents_map[player_index].add(opponent_index)
                    print()
            for player_index in range(player_count):
                assert len(player_opponents_map[player_index]) == 3 * round_count
            line_index += round_count
            print()


if __name__ == "__main__":
    main()
