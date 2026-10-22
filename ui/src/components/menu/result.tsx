import {
  columnSizingFeature,
  createColumnHelper,
  tableFeatures,
  useTable,
} from "@tanstack/react-table";
import {
  MDBContainer as Container,
  MDBTable as Table,
  MDBTableBody as TableBody,
  MDBTableHead as TableHead,
} from "mdbreact";
import { Fragment, useMemo, useState } from "react";
import type { ResultData } from "../../api";
import LoadingScreen from "../shared/loading-screen";

interface ResultProps {
  mode: "CHAR" | "DNA" | "PRO";
  status: "idle" | "pending" | "error" | "success";
  data?: ResultData;
}

type Row = { label: string; value: string | number | boolean };

const features = tableFeatures({ columnSizingFeature });
const helper = createColumnHelper<typeof features, Row>();
const columns = helper.columns([
  helper.accessor("label", { header: "", size: 180 }),
  helper.accessor("value", {
    header: "",
    cell: info => String(info.getValue()),
  }),
]);

function fixed(width: number) {
  return { width, minWidth: width };
}

function ResultTable({
  title,
  record,
}: {
  title?: string;
  record: Record<string, string | number | boolean>;
}) {
  const data = useMemo(
    () => Object.entries(record).map(([label, value]) => ({ label, value })),
    [record],
  );
  const table = useTable({ features, columns, data });

  return (
    <Table responsive bordered small hover className="text-nowrap">
      {title && (
        <TableHead>
          <tr>
            {table.getFlatHeaders().map(header => (
              <th
                key={header.id}
                scope="col"
                style={
                  header.column.id === "label"
                    ? fixed(header.column.getSize())
                    : undefined
                }
              >
                {header.column.id === "value" ? title : null}
              </th>
            ))}
          </tr>
        </TableHead>
      )}
      <TableBody>
        {table.getRowModel().rows.map(row => (
          <tr key={row.id}>
            {row.getAllCells().map(cell =>
              cell.column.id === "label" ? (
                <th key={cell.id} scope="row" style={fixed(cell.column.getSize())}>
                  <table.FlexRender cell={cell} />
                </th>
              ) : (
                <td key={cell.id}>
                  <table.FlexRender cell={cell} />
                </td>
              ),
            )}
          </tr>
        ))}
      </TableBody>
    </Table>
  );
}

function isMessage(data: ResultData): data is { data: string } {
  return typeof data === "object" && typeof data.data === "string";
}

function Result({ mode, status, data }: ResultProps) {
  const [viewLimit, setViewLimit] = useState(10);

  function handleViewMore() {
    setViewLimit(viewLimit + 10);
  }

  if (status === "idle") return null;
  if (status === "pending") return <LoadingScreen />;

  if (status === "success" && data) {
    const results = data;
    if (isMessage(results) || Object.keys(results).length === 0) {
      return (
        <Container id="results" className="mb-5">
          <h2 className="h2-responsive mx-md-2 my-4 text-center">
            {isMessage(results) ? results.data : null}
          </h2>
        </Container>
      );
    }

    if (mode === "CHAR")
      return (
        <Container id="results" className="mb-5">
          <ResultTable
            record={Object.fromEntries(
              Object.entries(results).map(([key, row]) => [key, row["1"]]),
            )}
          />
        </Container>
      );

    return (
      <Container id="results" className="mb-5">
        <h2 className="h2-responsive mx-md-2 my-4">
          {Object.keys(results).length} result
          {Object.keys(results).length !== 1 && "s"}
        </h2>
        {Object.keys(results).map(
          (item, i) =>
            i < viewLimit && (
              <Fragment key={item}>
                <ResultTable title={`Primer ${i + 1}`} record={results[item]} />
                {i + 1 === viewLimit && i + 1 < Object.keys(results).length && (
                  <div className="mt-4 mb-5 text-center">
                    <button
                      className="btn btn-info"
                      onClick={handleViewMore}
                      type="button"
                    >
                      View more...
                    </button>
                  </div>
                )}
              </Fragment>
            ),
        )}
      </Container>
    );
  }

  return (
    <Container id="results" className="mb-5 text-center">
      <p>
        Oops! Something went wrong on the server. Please try again later, or{" "}
        <a
          href="https://github.com/kvdomingo/primerdriver/issues"
          target="_blank"
          rel="noopener noreferrer"
        >
          report an issue
        </a>
        .
      </p>
    </Container>
  );
}

export default Result;
