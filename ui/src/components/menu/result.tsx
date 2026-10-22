import {
  MDBContainer as Container,
  MDBTable as Table,
  MDBTableBody as TableBody,
  MDBTableHead as TableHead,
} from "mdbreact";
import { Fragment, useState } from "react";
import type { ResultData } from "../../api";
import LoadingScreen from "../shared/loading-screen";

interface ResultProps {
  mode: "CHAR" | "DNA" | "PRO";
  status: "idle" | "pending" | "error" | "success";
  data?: ResultData;
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
          <Table responsive bordered small hover className="text-nowrap">
            <TableBody>
              {Object.keys(results).map(key => (
                <tr key={key}>
                  <th scope="row">{key}</th>
                  <td>{results[key]["1"]}</td>
                </tr>
              ))}
            </TableBody>
          </Table>
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
                <Table responsive className="text-nowrap" small bordered hover>
                  <TableHead>
                    <tr>
                      <th scope="col" />
                      <th scope="col">{`Primer ${i + 1}`}</th>
                    </tr>
                  </TableHead>
                  <TableBody>
                    {Object.keys(results[(i + 1).toString()]).map(key => (
                      <tr key={key}>
                        <th scope="row">{key}</th>
                        <td>{results[(i + 1).toString()][key]}</td>
                      </tr>
                    ))}
                  </TableBody>
                </Table>
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
